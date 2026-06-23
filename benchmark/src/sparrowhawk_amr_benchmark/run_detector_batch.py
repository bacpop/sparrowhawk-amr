from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

try:
    from .common import ensure_detector_binary, ensure_dir, read_csv, run_and_time, write_csv
except ImportError:
    from common import ensure_detector_binary, ensure_dir, read_csv, run_and_time, write_csv


def run_build_index(
    detector_bin: Path,
    db_root: Path,
    index_path: Path,
    k: int,
    alphabet: str,
    min_exact_gene_kmers: int,
    min_hierarchy_unit_kmers: int,
) -> None:
    ensure_dir(index_path.parent)
    result = run_and_time(
        [
            str(detector_bin),
            "index",
            "build",
            "--db-dir",
            str(db_root),
            "--out",
            str(index_path),
            "--alphabet",
            alphabet,
            "--k",
            str(k),
            "--min-exact-gene-kmers",
            str(min_exact_gene_kmers),
            "--min-hierarchy-unit-kmers",
            str(min_hierarchy_unit_kmers),
        ]
    )
    if result["returncode"] != 0:
        raise RuntimeError(
            result["stderr"] or result["stdout"] or f"index build failed for {alphabet} k={k}"
        )


def run_report_map(detector_bin: Path, index_path: Path, report_map_path: Path) -> None:
    ensure_dir(report_map_path.parent)
    result = run_and_time(
        [
            str(detector_bin),
            "index",
            "report-map",
            "--index",
            str(index_path),
            "--out",
            str(report_map_path),
        ]
    )
    if result["returncode"] != 0:
        raise RuntimeError(
            result["stderr"] or result["stdout"] or f"report-map failed for {index_path}"
        )


def run_one(
    row: dict[str, str],
    detector_bin: Path,
    index_path: Path,
    result_dir: Path,
    min_gene_fraction: float,
    min_gene_group_fraction: float,
    mode: str,
    gene_call_root: Path,
) -> dict[str, str]:
    assembly_id = row["assembly_id"]
    fasta_path = Path(row["local_fasta_path"])
    out_json = result_dir / f"{assembly_id}.json"
    if row.get("fetch_status") not in {"downloaded", "cached"} or not fasta_path.exists():
        return {
            **row,
            "mode": mode,
            "detector_status": "skipped_missing_fasta",
            "detector_json": "",
            "elapsed_seconds": "0",
        }
    if out_json.exists():
        return {
            **row,
            "mode": mode,
            "detector_status": "cached",
            "detector_json": str(out_json),
            "elapsed_seconds": "0",
        }

    if mode == "direct":
        command = [
            str(detector_bin),
            "detect",
            "direct",
            "--index",
            str(index_path),
            "--fasta",
            str(fasta_path),
            "--sample-name",
            assembly_id,
        ]
    elif mode == "cds":
        command = [
            str(detector_bin),
            "detect",
            "cds",
            "--index",
            str(index_path),
            "--assembly",
            str(fasta_path),
            "--out-dir",
            str(ensure_dir(gene_call_root / mode / assembly_id)),
            "--sample-name",
            assembly_id,
        ]
    elif mode == "protein_cds":
        command = [
            str(detector_bin),
            "detect",
            "cds",
            "--protein",
            "--index",
            str(index_path),
            "--assembly",
            str(fasta_path),
            "--out-dir",
            str(ensure_dir(gene_call_root / mode / assembly_id)),
            "--sample-name",
            assembly_id,
        ]
    else:
        raise ValueError(f"unknown detector mode: {mode}")

    command.extend(
        [
            "--min-gene-fraction",
            str(min_gene_fraction),
            "--min-gene-group-fraction",
            str(min_gene_group_fraction),
        ]
    )
    result = run_and_time(command)
    if result["returncode"] == 0:
        out_json.write_text(result["stdout"])
    return {
        **row,
        "mode": mode,
        "detector_status": "ok" if result["returncode"] == 0 and out_json.exists() else "failed",
        "detector_json": str(out_json) if out_json.exists() else "",
        "elapsed_seconds": str(result["elapsed_seconds"]),
    }


def parse_csv_ints(raw: str) -> list[int]:
    return [int(part.strip()) for part in raw.split(",") if part.strip()]


def parse_csv_floats(raw: str) -> list[float]:
    return [float(part.strip()) for part in raw.split(",") if part.strip()]


def parse_csv_modes(raw: str) -> list[str]:
    modes = []
    for part in raw.split(","):
        mode = part.strip().replace("-", "_")
        if not mode:
            continue
        if mode not in {"direct", "cds", "protein_cds"}:
            raise ValueError(f"unknown mode: {part}")
        modes.append(mode)
    return modes


def format_fraction_label(value: float) -> str:
    scaled = int(round(value * 1000))
    whole = scaled // 1000
    frac = scaled % 1000
    return f"{whole}p{frac:03d}"


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Rust detector on benchmark cohort")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--detector-manifest", type=Path, required=True)
    parser.add_argument("--detector-bin", type=Path)
    parser.add_argument("--db-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--modes", type=str, default="direct")
    parser.add_argument("--ks", type=str, default="15,17,21,31")
    parser.add_argument("--protein-ks", type=str, default="5")
    parser.add_argument("--min-gene-fractions", type=str, default="0.02,0.05,0.10")
    parser.add_argument(
        "--min-report-unit-fractions",
        "--min-gene-group-fractions",
        "--min-family-fractions",
        dest="min_gene_group_fractions",
        type=str,
        default="0.10,0.20,0.30",
    )
    parser.add_argument("--min-exact-gene-kmers", type=int, default=20)
    parser.add_argument("--min-hierarchy-unit-kmers", type=int, default=20)
    parser.add_argument("--protein-min-exact-gene-kmers", type=int, default=5)
    parser.add_argument("--protein-min-hierarchy-unit-kmers", type=int, default=5)
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()

    detector_bin = args.detector_bin if args.detector_bin else ensure_detector_binary(args.detector_manifest)
    if not detector_bin.exists():
        raise FileNotFoundError(f"detector binary not found: {detector_bin}")
    manifest_rows = read_csv(args.manifest)
    modes = parse_csv_modes(args.modes)
    ks = parse_csv_ints(args.ks)
    protein_ks = parse_csv_ints(args.protein_ks)
    gene_fractions = parse_csv_floats(args.min_gene_fractions)
    gene_group_fractions = parse_csv_floats(args.min_gene_group_fractions)

    for mode in modes:
        mode_ks = protein_ks if mode == "protein_cds" else ks
        alphabet = "protein" if mode == "protein_cds" else "dna"
        min_exact = (
            args.protein_min_exact_gene_kmers
            if mode == "protein_cds"
            else args.min_exact_gene_kmers
        )
        min_hierarchy = (
            args.protein_min_hierarchy_unit_kmers
            if mode == "protein_cds"
            else args.min_hierarchy_unit_kmers
        )
        for k in mode_ks:
            index_path = ensure_dir(args.out_dir / "indexes") / f"{alphabet}_k{k}.amridx"
            run_build_index(detector_bin, args.db_root, index_path, k, alphabet, min_exact, min_hierarchy)
            run_report_map(
                detector_bin,
                index_path,
                ensure_dir(args.out_dir / "report_maps") / f"{alphabet}_k{k}.tsv",
            )
            for min_gene in gene_fractions:
                gene_label = format_fraction_label(min_gene)
                for min_family in gene_group_fractions:
                    gene_group_label = format_fraction_label(min_family)
                    result_dir = ensure_dir(
                        args.out_dir / mode / f"k_{k}" / f"gf_{gene_label}__ff_{gene_group_label}"
                    )
                    results = []
                    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                        futures = [
                            pool.submit(
                                run_one,
                                row,
                                detector_bin,
                                index_path,
                                result_dir,
                                min_gene,
                                min_family,
                                mode,
                                ensure_dir(args.out_dir / "gene_calls"),
                            )
                            for row in manifest_rows
                        ]
                        for future in as_completed(futures):
                            results.append(future.result())
                    ordered = sorted(results, key=lambda row: (row["species"], row["assembly_id"]))
                    write_csv(
                        args.out_dir
                        / mode
                        / f"k_{k}"
                        / f"gf_{gene_label}__ff_{gene_group_label}_status.csv",
                        list(ordered[0].keys()) if ordered else [],
                        ordered,
                    )


if __name__ == "__main__":
    main()
