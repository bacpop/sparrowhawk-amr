#!/usr/bin/env bash

# Shared configuration for the AMRFinderPlus-native benchmark evaluation.
# Override any of these by exporting the variable before running the scripts.

set -euo pipefail

COMMAND_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARK_DIR="$(cd "$COMMAND_DIR/.." && pwd)"

DETECTOR_BIN="${DETECTOR_BIN:-$(cd "$BENCHMARK_DIR/.." && pwd)/target/release/sparrowhawk-amr}"
RUN="${RUN:-/mnt/datapool/vrbouza/projects/assembler_development/tests/2026_05_19_newevalamr}"
AMRFINDER_DB_ROOT="${AMRFINDER_DB_ROOT:-$RUN/amrfinderplus_native_db}"
DB="${DB:-$AMRFINDER_DB_ROOT/latest}"
AMRFINDER_DB_URL="${AMRFINDER_DB_URL:-https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest}"

DATA_DIR="${DATA_DIR:-$BENCHMARK_DIR/data}"
AMR_RECORDS="${AMR_RECORDS:-$DATA_DIR/amr_records.csv}"
SELECTED_MANIFEST="${SELECTED_MANIFEST:-$RUN/manifests/selected_assemblies.csv}"
SELECTION_SUMMARY="${SELECTION_SUMMARY:-$RUN/manifests/selection_summary.json}"
FETCHED_MANIFEST="${FETCHED_MANIFEST:-$RUN/manifests/fetched_assemblies.csv}"
FASTA_DIR="${FASTA_DIR:-$RUN/fasta}"

AMRFINDER_BIN="${AMRFINDER_BIN:-amrfinder}"
AMRFINDER_UPDATE_BIN="${AMRFINDER_UPDATE_BIN:-}"
JOBS="${JOBS:-80}"
TARGET_SIZE="${TARGET_SIZE:-2000}"

# NCBI Entrez is used to resolve assembly download URLs. Keep this lower than
# detector parallelism to avoid API throttling on large cohorts.
FETCH_JOBS="${FETCH_JOBS:-2}"
FETCH_RETRY_ROUNDS="${FETCH_RETRY_ROUNDS:-3}"
FETCH_RETRY_JOBS="${FETCH_RETRY_JOBS:-1}"
FETCH_RETRY_SLEEP="${FETCH_RETRY_SLEEP:-30}"

# Detector defaults. Keep this to one k value; no sweep.
MODES="${MODES:-direct}"
KS="${KS:-31}"
PROTEIN_KS="${PROTEIN_KS:-5}"
MIN_GENE_FRACTIONS="${MIN_GENE_FRACTIONS:-0.10}"
MIN_REPORT_UNIT_FRACTIONS="${MIN_REPORT_UNIT_FRACTIONS:-0.10}"
MIN_EXACT_GENE_KMERS="${MIN_EXACT_GENE_KMERS:-20}"
MIN_HIERARCHY_UNIT_KMERS="${MIN_HIERARCHY_UNIT_KMERS:-20}"
PROTEIN_MIN_EXACT_GENE_KMERS="${PROTEIN_MIN_EXACT_GENE_KMERS:-5}"
PROTEIN_MIN_HIERARCHY_UNIT_KMERS="${PROTEIN_MIN_HIERARCHY_UNIT_KMERS:-5}"

DETECTOR_OUT="${DETECTOR_OUT:-$RUN/detector}"
NATIVE_OUT="${NATIVE_OUT:-$RUN/native_amrfinder_plus}"
NATIVE_EFFECTIVE_OUT="${NATIVE_EFFECTIVE_OUT:-}"
NATIVE_AMRFINDER_VERSION="${NATIVE_AMRFINDER_VERSION:-}"
NATIVE_DB_VERSION="${NATIVE_DB_VERSION:-}"
COMPARISON_OUT="${COMPARISON_OUT:-$RUN/comparisons}"
FAILURE_OUT="${FAILURE_OUT:-$RUN/failure_analysis}"
INCLUDE_TYPES="${INCLUDE_TYPES:-AMR,STRESS,VIRULENCE}"

safe_path_component() {
  local value="$1"
  value="$(printf '%s' "$value" | tr -cs '[:alnum:]._:-' '_')"
  value="${value##_}"
  value="${value%%_}"
  printf '%s\n' "${value:-unknown}"
}

first_existing_line() {
  local path
  for path in "$@"; do
    if [[ -s "$path" ]]; then
      head -n 1 "$path"
      return
    fi
  done
  printf '%s\n' "unknown"
}

compute_native_effective_out() {
  local requested_out="$NATIVE_EFFECTIVE_OUT"
  local amrfinder_version
  local db_version
  amrfinder_version="$("$AMRFINDER_BIN" --version 2>&1 | head -n 1 || true)"
  db_version="$(first_existing_line "$DB/version.txt" "$DB/database_format_version.txt")"
  NATIVE_AMRFINDER_VERSION="$amrfinder_version"
  NATIVE_DB_VERSION="$db_version"
  if [[ -z "$requested_out" ]]; then
    NATIVE_EFFECTIVE_OUT="$NATIVE_OUT/amrfinder_$(safe_path_component "$amrfinder_version")__db_$(safe_path_component "$db_version")"
  fi
  export NATIVE_EFFECTIVE_OUT NATIVE_AMRFINDER_VERSION NATIVE_DB_VERSION
}

# Selection defaults. Set FULL_SET=1 to evaluate every assembly in AMR_RECORDS.
FULL_SET="${FULL_SET:-0}"
ESKAPEE_FLOOR="${ESKAPEE_FLOOR:-300}"
CLASS_FLOOR="${CLASS_FLOOR:-300}"
ESKAPEE_SPECIES="${ESKAPEE_SPECIES:-}"
ANTIBIOTIC_CLASSES="${ANTIBIOTIC_CLASSES:-}"
