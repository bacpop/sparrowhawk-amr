#!/usr/bin/env bash

# Compare Sparrowhawk-AMR against native AMRFinderPlus --plus, then generate
# failure-analysis files and index statistics.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/00_config.sh"

compute_native_effective_out
NATIVE_STATUS="$NATIVE_EFFECTIVE_OUT/amrfinder_status.csv"
HIERARCHY="$DB/ReferenceGeneHierarchy.txt"
REPORT_MAP_ROOT="$DETECTOR_OUT/report_maps"
UNIT_STATS_ROOT="$DETECTOR_OUT/unit_stats"

if [[ ! -f "$HIERARCHY" ]]; then
  echo "Missing required AMRFinderPlus hierarchy file: $HIERARCHY" >&2
  echo "Run 01_prepare_2000_dataset.sh to download the complete AMRFinderPlus DB metadata." >&2
  exit 1
fi

if [[ ! -f "$NATIVE_STATUS" ]]; then
  echo "Missing native AMRFinderPlus status CSV: $NATIVE_STATUS" >&2
  echo "Native AMRFinderPlus cache root: $NATIVE_OUT" >&2
  echo "Native AMRFinderPlus effective output: $NATIVE_EFFECTIVE_OUT" >&2
  echo "Run 02_run_detectors_and_native_plus.sh first." >&2
  exit 1
fi

if [[ ! -d "$REPORT_MAP_ROOT" ]]; then
  echo "Missing detector report maps: $REPORT_MAP_ROOT" >&2
  echo "Run 02_run_detectors_and_native_plus.sh first." >&2
  exit 1
fi

mkdir -p "$COMPARISON_OUT" "$FAILURE_OUT"

STATUS_ARGS=()
while IFS= read -r status_csv; do
  STATUS_ARGS+=(--status-csv "$status_csv")
done < <(find "$DETECTOR_OUT" -type f -name '*_status.csv' | sort)

if [[ ${#STATUS_ARGS[@]} -eq 0 ]]; then
  echo "No detector status CSV files found under: $DETECTOR_OUT" >&2
  echo "Run 02_run_detectors_and_native_plus.sh first." >&2
  exit 1
fi

echo "Comparing detector calls to native AMRFinderPlus --plus..."
uv --directory "$BENCHMARK_DIR" run amr-compare-amrfinder-batch \
  --amrfinder-status "$NATIVE_STATUS" \
  --detector-root "$DETECTOR_OUT" \
  --report-map-root "$REPORT_MAP_ROOT" \
  --unit-stats-root "$UNIT_STATS_ROOT" \
  --detector-bin "$DETECTOR_BIN" \
  --detector-root "$DETECTOR_OUT" \
  --db-dir "$DB" \
  --hierarchy "$HIERARCHY" \
  --include-types "$INCLUDE_TYPES" \
  --out-dir "$COMPARISON_OUT" \
  "${STATUS_ARGS[@]}"

echo "Writing summary report..."
uv --directory "$BENCHMARK_DIR" run amr-report-results \
  --selected-manifest "$SELECTED_MANIFEST" \
  --aggregate-metrics "$COMPARISON_OUT/aggregate_metrics.csv" \
  --out-md "$COMPARISON_OUT/summary.md"

echo "Analyzing failures..."
uv --directory "$BENCHMARK_DIR" run amr-analyze-failures \
  --comparison-dir "$COMPARISON_OUT" \
  --report-map-root "$REPORT_MAP_ROOT" \
  --unit-stats-root "$UNIT_STATS_ROOT" \
  --detector-bin "$DETECTOR_BIN" \
  --detector-root "$DETECTOR_OUT" \
  --db-dir "$DB" \
  --hierarchy "$HIERARCHY" \
  --include-types "$INCLUDE_TYPES" \
  --out-dir "$FAILURE_OUT"

echo "Creating plots..."
uv --directory "$BENCHMARK_DIR" run amr-plot-results \
  --aggregate-metrics "$COMPARISON_OUT/aggregate_metrics.csv" \
  --selected-manifest "$SELECTED_MANIFEST" \
  --species-metrics "$COMPARISON_OUT/species_metrics.csv" \
  --class-metrics "$COMPARISON_OUT/class_metrics.csv" \
  --species-class-metrics "$COMPARISON_OUT/species_class_metrics.csv" \
  --type-metrics "$COMPARISON_OUT/type_metrics.csv" \
  --out-dir "$COMPARISON_OUT/plots"

echo "Collecting index stats..."
for index in "$DETECTOR_OUT"/indexes/*.amridx; do
  [[ -e "$index" ]] || continue
  name="$(basename "$index" .amridx)"
  "$DETECTOR_BIN" index stats \
    --index "$index" \
    > "$RUN/index_stats_${name}.txt"
  "$DETECTOR_BIN" index unit-stats \
    --index "$index" \
    --db-dir "$DB" \
    --out "$RUN/index_unit_stats_${name}.tsv"
  du -h "$index" > "$RUN/index_size_${name}.txt"
done

echo "Main outputs:"
echo "  Aggregate metrics:      $COMPARISON_OUT/aggregate_metrics.csv"
echo "  Species metrics:        $COMPARISON_OUT/species_metrics.csv"
echo "  Class metrics:          $COMPARISON_OUT/class_metrics.csv"
echo "  Subclass metrics:       $COMPARISON_OUT/subclass_metrics.csv"
echo "  Type metrics:           $COMPARISON_OUT/type_metrics.csv"
echo "  Species-class metrics:  $COMPARISON_OUT/species_class_metrics.csv"
echo "  Summary:                $COMPARISON_OUT/summary.md"
echo "  Plots:                  $COMPARISON_OUT/plots"
echo "  Failure analysis:  $FAILURE_OUT"
echo "  Index stats:       $RUN/index_stats_*.txt"
echo "  Index sizes:       $RUN/index_size_*.txt"
