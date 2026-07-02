# sparrowhawk-amr benchmark

Python/uv orchestration for comparing `sparrowhawk-amr` detector calls against native AMRFinderPlus `--plus` results. The Rust crate remains the detector and index builder; cohort selection, FASTA fetching, native AMRFinderPlus runs, comparison, reporting, and failure analysis live here.

## Basic workflow

Build `sparrowhawk-amr` outside the benchmark, set paths in the environment, then run the command scripts from `benchmark/commands`:

```sh
# From the sparrowhawk-amr checkout, once per binary update:
cargo build --release

export RUN=/path/to/evaluation_outputs
# AMR_RECORDS defaults to benchmark/data/amr_records.csv
# export AMR_RECORDS=/path/to/input_amr_records.csv
export AMRFINDER_BIN=/path/to/amrfinder
export DETECTOR_BIN=/path/to/sparrowhawk-amr/target/release/sparrowhawk-amr

# Optional: evaluate every assembly in AMR_RECORDS instead of selecting a cohort.
# export FULL_SET=1

cd /path/to/sparrowhawk-amr/benchmark
./commands/01_prepare_dataset.sh
./commands/02_run_detectors_and_native_plus.sh
./commands/03_compare_analyze_and_index_stats.sh
```

The input CSV is bundled by default at `benchmark/data/amr_records.csv`, and `AMR_RECORDS` can still be overridden.

Important roots are independent: `FASTA_DIR`, `DETECTOR_OUT`, `NATIVE_OUT`, `COMPARISON_OUT`, and `FAILURE_OUT` can all point outside the code checkout.

The benchmark scripts derive their own benchmark directory from the command location; there is no need to export `BENCHMARK_ROOT`. They do not compile Rust code, so `DETECTOR_BIN` must point to an existing executable.

## Selection

By default, `amr-select-subset` does stratified selection with the existing richness score and balancing logic. It adds floors for ESKAPEE species and antibiotic classes, defaulting to 300 each, and lets the cohort grow beyond `TARGET_SIZE` if needed. If a requested floor is impossible from the input CSV, selection fails with a shortfall report.

Use `FULL_SET=1` or `amr-select-subset --full-set` to evaluate every assembly represented in the input CSV.

## Metrics

Public comparison outputs report only two metric levels:

- `exact`: strict exact allele/gene match.
- `report_unit`: the hierarchy-aware report unit used by the detector index.

Legacy `gene_group` and `family` names may still appear internally as compatibility aliases for detector JSON or old status filenames, but they are not emitted as public metric columns.
