from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

try:
    from .common import (
        decompress_gzip,
        download_to_path,
        ensure_dir,
        read_csv,
        resolve_assembly_download_urls,
        write_csv,
    )
except ImportError:
    from common import (
        decompress_gzip,
        download_to_path,
        ensure_dir,
        read_csv,
        resolve_assembly_download_urls,
        write_csv,
    )


def local_fasta_name(assembly_id: str) -> str:
    return f"{assembly_id}_genomic.fna"


def expected_fasta_path(row: dict[str, str], out_dir: Path) -> Path:
    explicit = row.get("local_fasta_path", "")
    if explicit:
        return Path(explicit)
    return out_dir / row["assembly_id"] / local_fasta_name(row["assembly_id"])


def fasta_validation_error(path: Path) -> str:
    if not path.exists():
        return "missing"
    if path.stat().st_size == 0:
        return "empty"
    with path.open("rt", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                return "" if stripped.startswith(">") else "not_fasta"
    return "empty"


def has_expected_fasta(row: dict[str, str], out_dir: Path) -> bool:
    return fasta_validation_error(expected_fasta_path(row, out_dir)) == ""


def inspect_existing_fastas(rows: list[dict[str, str]], out_dir: Path) -> dict[str, object]:
    existing = 0
    valid = 0
    invalid: list[dict[str, str]] = []
    for row in rows:
        path = expected_fasta_path(row, out_dir)
        error = fasta_validation_error(path)
        if path.exists():
            existing += 1
        if not error:
            valid += 1
        elif error != "missing":
            invalid.append(
                {
                    "assembly_id": row["assembly_id"],
                    "path": str(path),
                    "reason": error,
                }
            )
    return {"existing": existing, "valid": valid, "invalid": invalid}


def alternate_accessions(assembly_id: str) -> list[str]:
    accessions = [assembly_id]
    if assembly_id.startswith("GCA_"):
        accessions.append("GCF_" + assembly_id[4:])
    elif assembly_id.startswith("GCF_"):
        accessions.append("GCA_" + assembly_id[4:])
    return accessions


def fetch_one(row: dict[str, str], out_dir: Path) -> dict[str, str]:
    assembly_id = row["assembly_id"]
    assembly_dir = ensure_dir(out_dir / assembly_id)
    fasta_path = assembly_dir / local_fasta_name(assembly_id)
    gz_path = fasta_path.with_suffix(".fna.gz")
    expected_fasta = str(fasta_path)
    cache_error = fasta_validation_error(fasta_path)
    if not cache_error:
        return {
            **row,
            "local_fasta_path": expected_fasta,
            "download_url": row.get("download_url", ""),
            "fetch_status": "cached",
            "fetch_message": "",
        }
    if cache_error != "missing" and fasta_path.exists():
        fasta_path.unlink()

    attempts: list[str] = []
    for accession in alternate_accessions(assembly_id):
        try:
            urls = resolve_assembly_download_urls(accession)
        except Exception as exc:  # noqa: BLE001
            attempts.append(f"{accession}: resolve failed: {exc}")
            continue

        for url in urls:
            try:
                download_to_path(url, gz_path)
                decompress_gzip(gz_path, fasta_path)
                if gz_path.exists():
                    gz_path.unlink()
                validation_error = fasta_validation_error(fasta_path)
                if not validation_error:
                    return {
                        **row,
                        "local_fasta_path": expected_fasta,
                        "download_url": url,
                        "fetch_status": "downloaded",
                        "fetch_message": "",
                    }
                attempts.append(
                    f"{accession}: downloaded invalid FASTA at {expected_fasta}: {validation_error} from {url}"
                )
            except Exception as exc:  # noqa: BLE001
                attempts.append(f"{accession}: {url}: {exc}")
            finally:
                if gz_path.exists():
                    gz_path.unlink()

    return {
        **row,
        "local_fasta_path": expected_fasta,
        "download_url": "",
        "fetch_status": "failed",
        "fetch_message": f"expected_fasta={expected_fasta}; " + "; ".join(attempts),
    }


def cached_result(row: dict[str, str], out_dir: Path) -> dict[str, str]:
    fasta_path = expected_fasta_path(row, out_dir)
    return {
        **row,
        "local_fasta_path": str(fasta_path),
        "download_url": row.get("download_url", ""),
        "fetch_status": "cached",
        "fetch_message": "",
    }


def fetch_round(rows: list[dict[str, str]], out_dir: Path, jobs: int, label: str) -> list[dict[str, str]]:
    results = []
    total = len(rows)
    if total == 0:
        return results
    print(f"{label}: fetching {total} assembly FASTA file(s) with {max(1, jobs)} job(s).", file=sys.stderr)
    completed = 0
    downloaded = 0
    cached = 0
    failed = 0
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        futures = [pool.submit(fetch_one, row, out_dir) for row in rows]
        for future in as_completed(futures):
            row = future.result()
            results.append(row)
            completed += 1
            status = row.get("fetch_status", "")
            downloaded += status == "downloaded"
            cached += status == "cached"
            failed += status == "failed"
            if completed == total or completed % 25 == 0:
                print(
                    f"{label}: {completed}/{total} complete; downloaded={downloaded}; cached={cached}; failed={failed}",
                    file=sys.stderr,
                )
    return results


def merge_results(current: list[dict[str, str]], updates: list[dict[str, str]]) -> list[dict[str, str]]:
    by_assembly = {row["assembly_id"]: row for row in current}
    for row in updates:
        by_assembly[row["assembly_id"]] = row
    return list(by_assembly.values())


def missing_rows(rows: list[dict[str, str]], out_dir: Path) -> list[dict[str, str]]:
    return [row for row in rows if not has_expected_fasta(row, out_dir)]


def ordered_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(rows, key=lambda row: (row.get("species", ""), row["assembly_id"]))


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch assembly FASTAs for AMR benchmark cohort")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--retry-missing-rounds", type=int, default=3)
    parser.add_argument("--retry-jobs", type=int, default=1)
    parser.add_argument("--retry-sleep", type=float, default=30.0)
    parser.add_argument("--allow-missing", action="store_true")
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    rows = read_csv(args.manifest)
    if args.limit > 0:
        rows = rows[: args.limit]

    print(f"Checking existing FASTA files under {args.out_dir}...", file=sys.stderr)
    cache = inspect_existing_fastas(rows, args.out_dir)
    invalid = cache["invalid"]
    need_download = len(rows) - int(cache["valid"])
    print(
        f"Found {cache['existing']} expected FASTA path(s); {cache['valid']} valid; {len(invalid)} invalid/corrupt/empty.",
        file=sys.stderr,
    )
    print(f"Need to download {need_download} FASTA file(s).", file=sys.stderr)
    for item in invalid[:10]:
        print(
            f"Invalid cached FASTA: assembly_id={item['assembly_id']} reason={item['reason']} path={item['path']}",
            file=sys.stderr,
        )
    if len(invalid) > 10:
        print(f"... {len(invalid) - 10} additional invalid cached FASTA file(s) omitted.", file=sys.stderr)

    cached_rows = [cached_result(row, args.out_dir) for row in rows if has_expected_fasta(row, args.out_dir)]
    rows_to_fetch = [row for row in rows if not has_expected_fasta(row, args.out_dir)]
    results = cached_rows + fetch_round(rows_to_fetch, args.out_dir, args.jobs, "Initial fetch")
    for retry_round in range(1, args.retry_missing_rounds + 1):
        missing = missing_rows(results, args.out_dir)
        if not missing:
            break
        print(
            f"{len(missing)} expected FASTA files are still missing after fetch round {retry_round}; retrying with {args.retry_jobs} job(s).",
            file=sys.stderr,
        )
        if args.retry_sleep > 0:
            time.sleep(args.retry_sleep)
        results = merge_results(
            results,
            fetch_round(missing, args.out_dir, args.retry_jobs, f"Retry {retry_round}"),
        )

    final_missing = missing_rows(results, args.out_dir)
    ordered = ordered_rows(results)
    write_csv(args.out_csv, list(ordered[0].keys()) if ordered else [], ordered)

    missing_csv = args.out_csv.with_name(args.out_csv.stem + "_missing_fastas.csv")
    if final_missing:
        missing_ordered = ordered_rows(final_missing)
        write_csv(missing_csv, list(missing_ordered[0].keys()), missing_ordered)
        print(
            f"{len(final_missing)} expected FASTA files are still missing after all retries; details: {missing_csv}",
            file=sys.stderr,
        )
        if not args.allow_missing:
            raise SystemExit(1)
    elif missing_csv.exists():
        missing_csv.unlink()


if __name__ == "__main__":
    main()
