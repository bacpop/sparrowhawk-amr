from __future__ import annotations

import argparse
import collections
import csv
import json
import statistics
from pathlib import Path
from typing import Any

try:
    from .compare_amrfinder_batch import (
        gene_group_key,
        load_hierarchy,
        load_report_map_rows,
        report_map_lookup,
        reportable_rows,
        report_map_path,
    )
    from .common import read_csv, write_csv
except ImportError:
    from compare_amrfinder_batch import (
        gene_group_key,
        load_hierarchy,
        load_report_map_rows,
        report_map_lookup,
        reportable_rows,
        report_map_path,
    )
    from common import read_csv, write_csv


def split_set(raw: str) -> set[str]:
    return {part for part in raw.split(";") if part}


def load_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


## To get from AMRF+ results
def report_unit_for_native(row: dict[str, str], report_map: dict[str, str], hierarchy: dict[str, dict[str, str]]) -> str:
    accession = row.get("Closest reference accession", "")
    symbol = row.get("Element symbol", "")
    node = row.get("Hierarchy node", "")
    unit = report_map.get(accession) or report_map.get(symbol) or report_map.get(node)
    if unit:
        return unit
    if node:
        return f"hierarchy_node:{gene_group_key(node, hierarchy)}"
    return f"exact_gene:{symbol}"

def report_unit_for_detector(hit: dict[str, Any]) -> str:
    unit_id = hit.get("unit_id")
    if not unit_id:
        return ""
    if hit.get("unit_type"):
        return f"{hit['unit_type']}:{unit_id}"
    if hit.get("call_type") == "gene":
        return f"exact_gene:{unit_id}"
    if hit.get("call_type") in {"gene_group", "family"}:
        return f"hierarchy_node:{unit_id}"
    return ""


def parse_mode_file(path: Path) -> tuple[str, int]:
    # direct_k_31_fraction_gene_0.100_gene_group_0.100_assemblies.csv
    parts = path.name.split("_")
    if path.name.startswith("protein_cds_"):
        return "protein_cds", int(parts[3])
    return parts[0], int(parts[2])


def pct(values: list[float], threshold: float) -> float:
    if not values:
        return 0.0
    return sum(1 for value in values if value <= threshold) / len(values)


def median(values: list[float]) -> float:
    return statistics.median(values) if values else 0.0


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize native AMRFinderPlus comparison failures")
    parser.add_argument("--comparison-dir", type=Path, required=True)
    parser.add_argument("--report-map-root", type=Path, required=True)
    parser.add_argument("--hierarchy", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--include-types", default="AMR,STRESS,VIRULENCE")
    args = parser.parse_args()

    included_types = {value.strip().upper() for value in args.include_types.split(",") if value.strip()}
    hierarchy = load_hierarchy(args.hierarchy)
    per_mode_rows = []
    missed_unit_rows = []
    detector_only_rows = []

    for assembly_csv in sorted(args.comparison_dir.glob("*_assemblies.csv")):
        mode, k = parse_mode_file(assembly_csv)
        report_map = report_map_lookup(load_report_map_rows(report_map_path(args.report_map_root, mode, k)))
        rows = read_csv(assembly_csv)
        native_method_counts: collections.Counter[str] = collections.Counter()
        missed_units: collections.Counter[str] = collections.Counter()
        detector_units: collections.Counter[str] = collections.Counter()
        exact_resolved_by_report_unit = 0
        total_exact_fn = 0
        total_report_unit_fn = 0
        total_report_unit_fp = 0
        fp_fractions: list[float] = []
        fp_distinct: list[float] = []
        fp_diag: list[float] = []

        for row in rows:
            total_exact_fn += int(row["exact_fn"])
            total_report_unit_fn += int(row["report_unit_fn"])
            total_report_unit_fp += int(row["report_unit_fp"])
            exact_resolved_by_report_unit += max(0, int(row["exact_fn"]) - int(row["report_unit_fn"]))

            missed = split_set(row["baseline_only_report_unit"])
            detector_only = split_set(row["detector_only_report_unit"])
            missed_units.update(missed)
            detector_units.update(detector_only)

            native_rows = reportable_rows(Path(row["native_amrfinder_tsv"]), included_types)
            for native_row in native_rows:
                unit = report_unit_for_native(native_row, report_map, hierarchy)
                if unit in missed:
                    native_method_counts[native_row.get("Method", "")] += 1
                    missed_unit_rows.append(
                        {
                            "mode": mode,
                            "k": k,
                            "assembly_id": row["assembly_id"],
                            "species": row.get("species", ""),
                            "antibiotic_classes": row.get("classes", ""),
                            "report_unit": unit,
                            "element_symbol": native_row.get("Element symbol", ""),
                            "hierarchy_node": native_row.get("Hierarchy node", ""),
                            "method": native_row.get("Method", ""),
                            "coverage": native_row.get("% Coverage of reference", ""),
                            "identity": native_row.get("% Identity to reference", ""),
                            "closest_reference_accession": native_row.get("Closest reference accession", ""),
                        }
                    )

            detector_payload = load_json(Path(row["detector_json"]))
            for hit in detector_payload.get("hits", []):
                unit = report_unit_for_detector(hit)
                if not unit:
                    continue
                if unit in detector_only:
                    detector_only_rows.append(
                        {
                            "mode": mode,
                            "k": k,
                            "assembly_id": row["assembly_id"],
                            "species": row.get("species", ""),
                            "antibiotic_classes": row.get("classes", ""),
                            "report_unit": unit,
                            "call_type": hit.get("call_type", hit.get("unit_type", "")),
                            "unit_label": hit.get("unit_label", ""),
                            "query_id": hit.get("query_id", ""),
                            "call_fraction": hit.get("call_fraction", 0.0),
                            "first_pass_distinct": hit.get("first_pass_distinct", 0),
                            "first_pass_diagnostic_total": hit.get("first_pass_diagnostic_total", 0),
                            "class_name": hit.get("class_name", ""),
                            "subclass": hit.get("subclass", ""),
                        }
                    )
                    fp_fractions.append(float(hit.get("call_fraction", 0.0)))
                    fp_distinct.append(float(hit.get("first_pass_distinct", 0)))
                    fp_diag.append(float(hit.get("first_pass_diagnostic_total", 0)))

        per_mode_rows.append(
            {
                "mode": mode,
                "k": k,
                "assemblies": len(rows),
                "exact_fn": total_exact_fn,
                "report_unit_fn": total_report_unit_fn,
                "exact_fn_resolved_by_report_unit": exact_resolved_by_report_unit,
                "report_unit_fp": total_report_unit_fp,
                "top_missed_units": ";".join(f"{unit}:{count}" for unit, count in missed_units.most_common(10)),
                "missed_native_methods": ";".join(f"{method}:{count}" for method, count in native_method_counts.most_common()),
                "top_detector_only_units": ";".join(f"{unit}:{count}" for unit, count in detector_units.most_common(10)),
                "detector_only_hit_count": len(fp_fractions),
                "detector_only_call_fraction_median": median(fp_fractions),
                "detector_only_distinct_median": median(fp_distinct),
                "detector_only_diag_total_median": median(fp_diag),
                "detector_only_pct_distinct_le_2": pct(fp_distinct, 2.0),
                "detector_only_pct_diag_total_le_10": pct(fp_diag, 10.0),
            }
        )

    write_csv(args.out_dir / "failure_summary.csv", list(per_mode_rows[0].keys()) if per_mode_rows else [], per_mode_rows)
    write_csv(args.out_dir / "missed_report_units.csv", list(missed_unit_rows[0].keys()) if missed_unit_rows else [], missed_unit_rows)
    write_csv(args.out_dir / "detector_only_report_units.csv", list(detector_only_rows[0].keys()) if detector_only_rows else [], detector_only_rows)


if __name__ == "__main__":
    main()
