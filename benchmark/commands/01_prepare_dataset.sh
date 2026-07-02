#!/usr/bin/env bash

# Fetch the AMRFinderPlus DB, select 2000 assemblies, and download their FASTA files.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/00_config.sh"

mkdir -p "$RUN"

if [[ ! -x "$DETECTOR_BIN" ]]; then
  echo "Missing DETECTOR_BIN: $DETECTOR_BIN" >&2
  echo "Set DETECTOR_BIN=/path/to/sparrowhawk-amr, or build the release binary before running this benchmark." >&2
  exit 1
fi

find_amrfinder_update() {
  if [[ -n "$AMRFINDER_UPDATE_BIN" ]]; then
    printf '%s\n' "$AMRFINDER_UPDATE_BIN"
    return
  fi
  if command -v amrfinder_update >/dev/null 2>&1; then
    command -v amrfinder_update
    return
  fi
  local amrfinder_dir
  amrfinder_dir="$(dirname "$AMRFINDER_BIN")"
  if [[ -x "$amrfinder_dir/amrfinder_update" ]]; then
    printf '%s\n' "$amrfinder_dir/amrfinder_update"
    return
  fi
  return 1
}

AMRFINDER_UPDATE_BIN_RESOLVED="$(find_amrfinder_update || true)"
if [[ -z "$AMRFINDER_UPDATE_BIN_RESOLVED" ]]; then
  echo "Missing amrfinder_update; cannot download and index the native AMRFinderPlus database." >&2
  echo "Set AMRFINDER_UPDATE_BIN=/path/to/amrfinder_update or load the AMRFinderPlus environment." >&2
  exit 1
fi

echo "Downloading and indexing AMRFinderPlus database under $AMRFINDER_DB_ROOT ..."
mkdir -p "$AMRFINDER_DB_ROOT"
"$AMRFINDER_UPDATE_BIN_RESOLVED" -d "$AMRFINDER_DB_ROOT"

if [[ ! -d "$DB" ]]; then
  echo "AMRFinderPlus update finished, but DB directory is missing: $DB" >&2
  echo "If you override DB, point it to the prepared database directory, usually AMRFINDER_DB_ROOT/latest." >&2
  exit 1
fi

download_db_file_if_missing() {
  local name="$1"
  local dest="$DB/$name"
  local url="${AMRFINDER_DB_URL%/}/$name"
  if [[ -s "$dest" ]]; then
    return
  fi
  echo "Downloading missing AMRFinderPlus DB file: $name"
  if command -v curl >/dev/null 2>&1; then
    curl -L -f --retry 3 --retry-delay 5 --retry-connrefused -o "$dest" "$url"
  elif command -v wget >/dev/null 2>&1; then
    wget --tries=3 --waitretry=5 -O "$dest" "$url"
  else
    echo "Missing curl/wget; cannot download $url" >&2
    exit 1
  fi
  if [[ ! -s "$dest" ]]; then
    echo "Downloaded DB file is missing or empty: $dest" >&2
    exit 1
  fi
}

# amrfinder_update prepares the BLAST databases, but some AMRFinderPlus metadata
# files can be absent in the local installation layout. Pull any missing flat
# files from the official latest database directory and fail if required files
# are still not present.
AMRFINDER_DB_FILES=(
  AMR.LIB
  AMRProt-mutation.tsv
  AMRProt-suppress.tsv
  AMRProt-susceptible.fa
  AMRProt-susceptible.tsv
  AMRProt.fa
  AMR_CDS.fa
  AMR_DNA-Acinetobacter_baumannii.fa
  AMR_DNA-Acinetobacter_baumannii.tsv
  AMR_DNA-Bordetella_pertussis.fa
  AMR_DNA-Bordetella_pertussis.tsv
  AMR_DNA-Campylobacter.fa
  AMR_DNA-Campylobacter.tsv
  AMR_DNA-Clostridioides_difficile.fa
  AMR_DNA-Clostridioides_difficile.tsv
  AMR_DNA-Enterococcus_faecalis.fa
  AMR_DNA-Enterococcus_faecalis.tsv
  AMR_DNA-Enterococcus_faecium.fa
  AMR_DNA-Enterococcus_faecium.tsv
  AMR_DNA-Escherichia.fa
  AMR_DNA-Escherichia.tsv
  AMR_DNA-Klebsiella_oxytoca.fa
  AMR_DNA-Klebsiella_oxytoca.tsv
  AMR_DNA-Klebsiella_pneumoniae.fa
  AMR_DNA-Klebsiella_pneumoniae.tsv
  AMR_DNA-Neisseria_gonorrhoeae.fa
  AMR_DNA-Neisseria_gonorrhoeae.tsv
  AMR_DNA-Salmonella.fa
  AMR_DNA-Salmonella.tsv
  AMR_DNA-Staphylococcus_aureus.fa
  AMR_DNA-Staphylococcus_aureus.tsv
  AMR_DNA-Streptococcus_pneumoniae.fa
  AMR_DNA-Streptococcus_pneumoniae.tsv
  ReferenceGeneCatalog.txt
  ReferenceGeneHierarchy.txt
  amr_targets.fa
  changelog.txt
  changes.txt
  database_format_version.txt
  fam.tsv
  mapgenelist.txt
  taxgroup.tsv
  version.txt
)
for db_file in "${AMRFINDER_DB_FILES[@]}"; do
  download_db_file_if_missing "$db_file"
done

REQUIRED_DB_FILES=(
  AMRProt.fa
  AMR_CDS.fa
  ReferenceGeneCatalog.txt
  ReferenceGeneHierarchy.txt
  fam.tsv
  mapgenelist.txt
  taxgroup.tsv
  version.txt
)
for required in "${REQUIRED_DB_FILES[@]}"; do
  if [[ ! -s "$DB/$required" ]]; then
    echo "Required AMRFinderPlus database file is missing after download: $DB/$required" >&2
    exit 1
  fi
done

if [[ ! -f "$DB/AMRProt.fa.phr" ]]; then
  echo "AMRFinderPlus database is missing the prepared BLAST protein database: $DB/AMRProt.fa.phr" >&2
  echo "Try forcing the AMRFinderPlus update: \"$AMRFINDER_UPDATE_BIN_RESOLVED\" --force_update -d \"$AMRFINDER_DB_ROOT\"" >&2
  exit 1
fi

if [[ ! -f "$AMR_RECORDS" ]]; then
  echo "Missing AMR_RECORDS input: $AMR_RECORDS" >&2
  echo "Set AMR_RECORDS=/path/to/amr_records.csv or provide this file." >&2
  exit 1
fi

echo "Selecting assemblies..."
mkdir -p "$(dirname "$SELECTED_MANIFEST")"
SELECT_ARGS=(
  --csv "$AMR_RECORDS"
  --out-csv "$SELECTED_MANIFEST"
  --out-json "$SELECTION_SUMMARY"
  --target-size "$TARGET_SIZE"
  --eskapee-floor "$ESKAPEE_FLOOR"
  --class-floor "$CLASS_FLOOR"
)
if [[ "$FULL_SET" == "1" ]]; then
  SELECT_ARGS+=(--full-set)
fi
if [[ -n "$ESKAPEE_SPECIES" ]]; then
  SELECT_ARGS+=(--eskapee-species "$ESKAPEE_SPECIES")
fi
if [[ -n "$ANTIBIOTIC_CLASSES" ]]; then
  SELECT_ARGS+=(--antibiotic-classes "$ANTIBIOTIC_CLASSES")
fi
uv --directory "$BENCHMARK_DIR" run amr-select-subset "${SELECT_ARGS[@]}"

echo "Fetching selected assemblies..."
mkdir -p "$(dirname "$FETCHED_MANIFEST")"
uv --directory "$BENCHMARK_DIR" run amr-fetch-assemblies \
  --manifest "$SELECTED_MANIFEST" \
  --out-dir "$FASTA_DIR" \
  --out-csv "$FETCHED_MANIFEST" \
  --jobs "$FETCH_JOBS" \
  --retry-missing-rounds "$FETCH_RETRY_ROUNDS" \
  --retry-jobs "$FETCH_RETRY_JOBS" \
  --retry-sleep "$FETCH_RETRY_SLEEP"

echo "Prepared:"
echo "  DB: $DB"
echo "  Selected manifest: $SELECTED_MANIFEST"
echo "  Fetched manifest: $FETCHED_MANIFEST"
