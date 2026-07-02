#!/usr/bin/env bash

# Run Sparrowhawk-AMR and native AMRFinderPlus --plus over the fetched manifest.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/00_config.sh"

if [[ ! -f "$FETCHED_MANIFEST" ]]; then
  echo "Missing fetched manifest: $FETCHED_MANIFEST" >&2
  echo "Run 01_prepare_2000_dataset.sh first, or set FETCHED_MANIFEST." >&2
  exit 1
fi

if [[ ! -d "$DB" ]]; then
  echo "Missing AMRFinderPlus DB: $DB" >&2
  echo "Run 01_prepare_2000_dataset.sh first, or set DB." >&2
  exit 1
fi

compute_native_effective_out
mkdir -p "$DETECTOR_OUT" "$NATIVE_EFFECTIVE_OUT"

echo "Running Sparrowhawk-AMR detector..."
uv --directory "$BENCHMARK_DIR" run amr-run-detector-batch \
  --manifest "$FETCHED_MANIFEST" \
  --detector-bin "$DETECTOR_BIN" \
  --db-root "$DB" \
  --out-dir "$DETECTOR_OUT" \
  --modes "$MODES" \
  --ks "$KS" \
  --protein-ks "$PROTEIN_KS" \
  --min-gene-fractions "$MIN_GENE_FRACTIONS" \
  --min-report-unit-fractions "$MIN_REPORT_UNIT_FRACTIONS" \
  --min-exact-gene-kmers "$MIN_EXACT_GENE_KMERS" \
  --min-hierarchy-unit-kmers "$MIN_HIERARCHY_UNIT_KMERS" \
  --protein-min-exact-gene-kmers "$PROTEIN_MIN_EXACT_GENE_KMERS" \
  --protein-min-hierarchy-unit-kmers "$PROTEIN_MIN_HIERARCHY_UNIT_KMERS" \
  --jobs "$JOBS"

echo "Running native AMRFinderPlus with --plus..."
export RUN FETCHED_MANIFEST NATIVE_OUT NATIVE_EFFECTIVE_OUT NATIVE_AMRFINDER_VERSION NATIVE_DB_VERSION AMRFINDER_BIN DB JOBS FASTA_DIR
python3 - <<'PY'
import csv
import json
import os
import pathlib
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

manifest = pathlib.Path(os.environ["FETCHED_MANIFEST"])
outdir = pathlib.Path(os.environ["NATIVE_EFFECTIVE_OUT"])
outdir.mkdir(parents=True, exist_ok=True)

amrfinder = os.environ["AMRFINDER_BIN"]
db = os.environ["DB"]
jobs = int(os.environ.get("JOBS", "8"))

metadata_path = outdir / "native_cache_metadata.json"
metadata_path.write_text(
    json.dumps(
        {
            "native_out_root": os.environ["NATIVE_OUT"],
            "native_effective_out": str(outdir),
            "amrfinder_bin": amrfinder,
            "amrfinder_version": os.environ.get("NATIVE_AMRFINDER_VERSION", ""),
            "db": db,
            "db_version": os.environ.get("NATIVE_DB_VERSION", ""),
        },
        indent=2,
        sort_keys=True,
    )
    + "\n"
)

with manifest.open() as fh:
    rows = list(csv.DictReader(fh))

def get_required(row, *names):
    for name in names:
        value = row.get(name)
        if value:
            return value
    raise KeyError(f"Missing one of {names} in manifest columns {list(row)}")


def get_optional(row, *names):
    for name in names:
        value = row.get(name)
        if value:
            return value
    return ""


def expected_fasta_path(row, assembly_id):
    explicit = get_optional(row, "local_fasta_path", "fasta_path", "fasta")
    if explicit:
        return pathlib.Path(explicit)
    return pathlib.Path(os.environ["FASTA_DIR"]) / assembly_id / f"{assembly_id}_genomic.fna"


def native_tsv_path(assembly_id):
    return outdir / f"{assembly_id}.amrfinder.tsv"


def valid_amrfinder_tsv(path):
    if not path.exists() or path.stat().st_size == 0:
        return False
    with path.open(newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    required = {"Element symbol", "Type", "Subtype", "Class", "Method"}
    return required.issubset(set(header))


def inspect_native_cache(rows):
    existing = 0
    valid = 0
    invalid = []
    for row in rows:
        assembly_id = get_required(row, "assembly_id", "assembly_accession", "accession")
        tsv_path = native_tsv_path(assembly_id)
        if tsv_path.exists():
            existing += 1
            if valid_amrfinder_tsv(tsv_path):
                valid += 1
            else:
                invalid.append({"assembly_id": assembly_id, "path": str(tsv_path)})
    return existing, valid, invalid


def skip_missing_fasta(assembly_id, fasta_path, row, message):
    (outdir / f"{assembly_id}.missing_fasta.txt").write_text(message + "\n")
    return {
        "assembly_id": assembly_id,
        "tsv_path": "",
        "returncode": 1,
        "seconds": "0.000",
        "expected_fasta": str(fasta_path),
        "fetch_status": row.get("fetch_status", ""),
        "native_status": "missing_fasta",
        "message": message,
    }


def run_one(row):
    assembly_id = get_required(row, "assembly_id", "assembly_accession", "accession")
    fasta_path = expected_fasta_path(row, assembly_id)
    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        message = (
            f"Missing FASTA for {assembly_id}; expected {fasta_path}. "
            f"fetch_status={row.get('fetch_status', '')}; fetch_message={row.get('fetch_message', '')}"
        )
        return skip_missing_fasta(assembly_id, fasta_path, row, message)

    tsv_path = native_tsv_path(assembly_id)
    if valid_amrfinder_tsv(tsv_path):
        return {
            "assembly_id": assembly_id,
            "tsv_path": str(tsv_path),
            "returncode": 0,
            "seconds": "0.000",
            "expected_fasta": str(fasta_path),
            "fetch_status": row.get("fetch_status", ""),
            "native_status": "cached",
            "message": "",
        }

    t0 = time.time()
    cmd = [amrfinder, "--plus", "-n", str(fasta_path), "-d", db, "-o", str(tsv_path)]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0:
        (outdir / f"{assembly_id}.stdout.txt").write_text(proc.stdout)
        (outdir / f"{assembly_id}.stderr.txt").write_text(proc.stderr)
    return {
        "assembly_id": assembly_id,
        "tsv_path": str(tsv_path) if proc.returncode == 0 else "",
        "returncode": proc.returncode,
        "seconds": f"{time.time() - t0:.3f}",
        "expected_fasta": str(fasta_path),
        "fetch_status": row.get("fetch_status", ""),
        "native_status": "ran" if proc.returncode == 0 else "failed",
        "message": "" if proc.returncode == 0 else "AMRFinderPlus failed; see stdout/stderr files",
    }

existing, valid, invalid = inspect_native_cache(rows)
need_run = len(rows) - valid
print(f"Checking cached native AMRFinderPlus TSVs under {outdir}...")
print(f"Found {existing} existing native TSV(s); {valid} valid; {len(invalid)} invalid/empty.")
print(f"Need to run AMRFinderPlus for {need_run} assembly/assemblies.")
for item in invalid[:10]:
    print(f"Invalid native TSV: assembly_id={item['assembly_id']} path={item['path']}")
if len(invalid) > 10:
    print(f"... {len(invalid) - 10} additional invalid native TSV file(s) omitted.")

status = []
completed = 0
cached = 0
ran = 0
failed = 0
missing_fasta = 0
with ThreadPoolExecutor(max_workers=jobs) as executor:
    futures = [executor.submit(run_one, row) for row in rows]
    for future in as_completed(futures):
        result = future.result()
        status.append(result)
        completed += 1
        native_status = result.get("native_status", "")
        cached += native_status == "cached"
        ran += native_status == "ran"
        failed += native_status == "failed"
        missing_fasta += native_status == "missing_fasta"
        if completed == len(rows) or completed % 25 == 0:
            print(
                f"AMRFinderPlus: {completed}/{len(rows)} complete; cached={cached}; ran={ran}; failed={failed}; missing_fasta={missing_fasta}"
            )

status_path = outdir / "amrfinder_status.csv"
with status_path.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=["assembly_id", "tsv_path", "returncode", "seconds", "expected_fasta", "fetch_status", "native_status", "message"],
    )
    writer.writeheader()
    writer.writerows(sorted(status, key=lambda row: row["assembly_id"]))

print(f"Wrote {status_path}")
PY

echo "Detector output: $DETECTOR_OUT"
echo "Native AMRFinderPlus cache root: $NATIVE_OUT"
echo "Native AMRFinderPlus output: $NATIVE_EFFECTIVE_OUT"
