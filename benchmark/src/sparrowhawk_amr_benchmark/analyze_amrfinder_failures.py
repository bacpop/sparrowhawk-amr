from __future__ import annotations

import argparse
import collections
import csv
import json
import re
import statistics
import subprocess
from pathlib import Path
from typing import Any

try:
    from .compare_amrfinder_batch import (
        canonical_report_node,
        detector_report_units_for_matching,
        filter_detector_hits_by_type,
        load_hierarchy,
        load_report_map_rows,
        report_map_context,
        report_map_lookup,
        report_unit_differences,
        reportable_rows,
        report_map_path,
        unit_stats_path,
        load_unit_stats_rows,
        unit_stats_lookup,
    )
    from .common import read_csv, write_csv
except ImportError:
    from compare_amrfinder_batch import (
        canonical_report_node,
        detector_report_units_for_matching,
        filter_detector_hits_by_type,
        load_hierarchy,
        load_report_map_rows,
        report_map_context,
        report_map_lookup,
        report_unit_differences,
        reportable_rows,
        report_map_path,
        unit_stats_path,
        load_unit_stats_rows,
        unit_stats_lookup,
    )
    from common import read_csv, write_csv



def load_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)



def resolve_path(raw: str, base: Path) -> Path:
    path = Path(raw)
    if path.is_absolute():
        return path
    candidate = base / path
    if candidate.exists():
        return candidate
    return path


def detector_index_path(detector_root: Path, mode: str, k: int) -> Path:
    alphabet = "protein" if mode == "protein_cds" else "dna"
    return detector_root / "indexes" / f"{alphabet}_k{k}.amridx"


def evidence_cache_path(out_dir: Path, assembly_csv: Path, assembly_id: str) -> Path:
    safe_id = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in assembly_id)
    return out_dir / "truth_kmer_evidence" / assembly_csv.stem / f"{safe_id}.json"


def run_truth_evidence(
    detector_bin: Path,
    detector_root: Path,
    db_dir: Path,
    assembly_csv: Path,
    row: dict[str, str],
    mode: str,
    k: int,
    include_types: set[str],
    out_dir: Path,
) -> list[dict[str, object]]:
    if mode == "protein_cds":
        return []
    assembly_fasta = row.get("assembly_fasta", "")
    if not assembly_fasta:
        return []
    out_json = evidence_cache_path(out_dir, assembly_csv, row["assembly_id"])
    out_json.parent.mkdir(parents=True, exist_ok=True)
    if not out_json.exists():
        cmd = [
            str(detector_bin),
            "eval",
            "truth-kmer-evidence",
            "--index",
            str(detector_index_path(detector_root, mode, k)),
            "--assembly",
            str(resolve_path(assembly_fasta, assembly_csv.parent)),
            "--amrfinder-tsv",
            str(resolve_path(row["native_amrfinder_tsv"], assembly_csv.parent)),
            "--detector-json",
            str(resolve_path(row["detector_json"], assembly_csv.parent)),
            "--db-dir",
            str(db_dir),
            "--include-types",
            ",".join(sorted(value.lower() for value in include_types)),
            "--min-gene-fraction",
            row.get("min_gene_fraction") or "0.100",
            "--min-family-fraction",
            row.get("min_report_unit_fraction") or "0.100",
            "--out",
            str(out_json),
        ]
        proc = subprocess.run(cmd, text=True, capture_output=True)
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr or proc.stdout or f"truth-kmer-evidence failed for {row['assembly_id']}")
    with out_json.open() as handle:
        payload = json.load(handle)
    return payload.get("rows", [])


def evidence_key(row: dict[str, object]) -> tuple[str, str, str]:
    return (
        str(row.get("element_symbol", "")),
        str(row.get("hierarchy_node", "")),
        str(row.get("closest_reference_accession", "")),
    )


def native_evidence_key(row: dict[str, str]) -> tuple[str, str, str]:
    return (
        row.get("Element symbol", ""),
        row.get("Hierarchy node", ""),
        row.get("Closest reference accession", ""),
    )


def evidence_fields(evidence: dict[str, object] | None) -> dict[str, object]:
    if not evidence:
        return {
            "truth_supported_by_index": "",
            "best_index_unit": "",
            "best_index_unit_type": "",
            "best_index_unit_label": "",
            "best_diagnostic_total": "",
            "best_diagnostic_matched": "",
            "best_diagnostic_missing": "",
            "best_diagnostic_fraction": "",
            "exact_diagnostic_total": "",
            "exact_diagnostic_matched": "",
            "exact_diagnostic_fraction": "",
            "family_diagnostic_total": "",
            "family_diagnostic_matched": "",
            "family_diagnostic_fraction": "",
            "interval_length": "",
            "interval_distinct_kmers": "",
            "recall_failure_category": "",
        }
    keys = [
        "truth_supported_by_index",
        "best_index_unit",
        "best_index_unit_type",
        "best_index_unit_label",
        "best_diagnostic_total",
        "best_diagnostic_matched",
        "best_diagnostic_missing",
        "best_diagnostic_fraction",
        "exact_diagnostic_total",
        "exact_diagnostic_matched",
        "exact_diagnostic_fraction",
        "family_diagnostic_total",
        "family_diagnostic_matched",
        "family_diagnostic_fraction",
        "interval_length",
        "interval_distinct_kmers",
        "recall_failure_category",
    ]
    return {key: evidence.get(key, "") for key in keys}


def numeric(value: object) -> float:
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return 0.0

def report_unit_for_native(row: dict[str, str], report_map: dict[str, str]) -> str:
    node = row.get("Hierarchy node", "")
    if node:
        return canonical_report_node(node)
    accession = row.get("Closest reference accession", "")
    symbol = row.get("Element symbol", "")
    unit = report_map.get(accession) or report_map.get(symbol)
    return canonical_report_node(unit or symbol)


def report_units_for_detector(hit: dict[str, Any]) -> set[str]:
    return detector_report_units_for_matching(hit)


def hit_total(hit: dict[str, Any], stats: dict[str, str] | None = None) -> str:
    value = hit.get("first_pass_diagnostic_total", "")
    if value not in (None, ""):
        return str(value)
    if stats:
        return stats.get("diagnostic_kmers", "")
    return ""


def hit_matched(hit: dict[str, Any]) -> str:
    value = hit.get("first_pass_distinct", "")
    return "" if value is None else str(value)


def missing_count(total: object, matched: object) -> str:
    try:
        return str(max(0, int(float(str(total))) - int(float(str(matched)))))
    except (TypeError, ValueError):
        return ""


def fraction(matched: object, total: object) -> str:
    try:
        numerator = float(str(matched))
        denominator = float(str(total))
    except (TypeError, ValueError):
        return ""
    if denominator <= 0:
        return "0.0"
    return f"{numerator / denominator:.6f}"


MODE_FILE_RE = re.compile(
    r"^(?P<mode>direct|cds|protein_cds)_k_(?P<k>\d+)_(?P<tmode>fraction|absolute)"
    r"_gene_(?P<gene>.+?)_(?:report_unit|gene_group)_(?P<unit>.+?)_assemblies$"
)


def parse_mode_file(path: Path) -> tuple[str, int, str, str, str]:
    # e.g. direct_k_31_fraction_gene_0.100_report_unit_0.100_assemblies.csv
    match = MODE_FILE_RE.fullmatch(path.stem)
    if not match:
        raise ValueError(f"Unexpected comparison filename: {path.name}")
    return match["mode"], int(match["k"]), match["tmode"], match["gene"], match["unit"]


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
    parser.add_argument("--unit-stats-root", type=Path)
    parser.add_argument("--detector-bin", type=Path)
    parser.add_argument("--detector-root", type=Path)
    parser.add_argument("--db-dir", type=Path)
    parser.add_argument("--hierarchy", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--include-types", default="AMR,STRESS,VIRULENCE")
    args = parser.parse_args()

    included_types = {value.strip().upper() for value in args.include_types.split(",") if value.strip()}
    hierarchy = load_hierarchy(args.hierarchy)
    per_mode_rows = []
    missed_unit_rows = []
    detector_only_rows = []
    truth_evidence_rows = []
    per_unit: dict[tuple[str, int, str, str, str], collections.Counter[str]] = collections.defaultdict(collections.Counter)

    for assembly_csv in sorted(args.comparison_dir.glob("*_assemblies.csv")):
        mode, k, _threshold_mode, gene_thr, unit_thr = parse_mode_file(assembly_csv)
        report_map_rows = load_report_map_rows(report_map_path(args.report_map_root, mode, k))
        report_map = report_map_lookup(report_map_rows)
        report_context = report_map_context(report_map_rows, hierarchy, included_types)
        if args.unit_stats_root:
            unit_stats = unit_stats_lookup(load_unit_stats_rows(unit_stats_path(args.unit_stats_root, mode, k)))
        else:
            unit_stats = {}
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

            native_rows = reportable_rows(Path(row["native_amrfinder_tsv"]), included_types)
            detector_payload = filter_detector_hits_by_type(load_json(Path(row["detector_json"])), included_types)
            evidence_by_key: dict[tuple[str, str, str], dict[str, object]] = {}
            if args.detector_bin and args.detector_root and args.db_dir:
                for evidence in run_truth_evidence(
                    args.detector_bin,
                    args.detector_root,
                    args.db_dir,
                    assembly_csv,
                    row,
                    mode,
                    k,
                    included_types,
                    args.out_dir,
                ):
                    evidence_row = {
                        "mode": mode,
                        "k": k,
                        "min_gene_threshold": gene_thr,
                        "min_report_unit_threshold": unit_thr,
                        "assembly_id": row["assembly_id"],
                        "species": row.get("species", ""),
                        "antibiotic_classes": row.get("classes", ""),
                        **evidence,
                    }
                    truth_evidence_rows.append(evidence_row)
                    evidence_by_key.setdefault(evidence_key(evidence), evidence)
            truth_units = {report_unit_for_native(native_row, report_map) for native_row in native_rows}
            detector_report_units = {
                unit
                for hit in detector_payload.get("hits", [])
                for unit in report_units_for_detector(hit)
            }
            missed, detector_only = report_unit_differences(
                detector_report_units,
                truth_units,
                report_context,
            )
            missed_units.update(missed)
            detector_units.update(detector_only)
            for unit in truth_units - missed:
                per_unit[(mode, k, gene_thr, unit_thr, unit)]["tp"] += 1
            for unit in missed:
                per_unit[(mode, k, gene_thr, unit_thr, unit)]["fn"] += 1
            for unit in detector_only:
                per_unit[(mode, k, gene_thr, unit_thr, unit)]["fp"] += 1
            for native_row in native_rows:
                unit = report_unit_for_native(native_row, report_map)
                if unit in missed:
                    native_method_counts[native_row.get("Method", "")] += 1
                    stats = unit_stats.get(unit, {})
                    evidence = evidence_by_key.get(native_evidence_key(native_row))
                    missed_unit_rows.append(
                        {
                            "mode": mode,
                            "k": k,
                            "min_gene_threshold": gene_thr,
                            "min_report_unit_threshold": unit_thr,
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
                            "unit_type": stats.get("unit_type", ""),
                            "diagnostic_total": stats.get("diagnostic_kmers", ""),
                            "diagnostic_matched": evidence_fields(evidence).get("best_diagnostic_matched", ""),
                            "diagnostic_missing": evidence_fields(evidence).get("best_diagnostic_missing", ""),
                            "diagnostic_fraction": evidence_fields(evidence).get("best_diagnostic_fraction", ""),
                            "member_genes": stats.get("member_genes", ""),
                            **evidence_fields(evidence),
                        }
                    )

            for hit in detector_payload.get("hits", []):
                for unit in report_units_for_detector(hit):
                    if unit not in detector_only:
                        continue
                    stats = unit_stats.get(unit, {})
                    total = hit_total(hit, stats)
                    matched = hit_matched(hit)
                    detector_only_rows.append(
                        {
                            "mode": mode,
                            "k": k,
                            "min_gene_threshold": gene_thr,
                            "min_report_unit_threshold": unit_thr,
                            "assembly_id": row["assembly_id"],
                            "species": row.get("species", ""),
                            "antibiotic_classes": row.get("classes", ""),
                            "report_unit": unit,
                            "call_type": hit.get("call_type", hit.get("unit_type", "")),
                            "unit_type": hit.get("unit_type", stats.get("unit_type", "")),
                            "unit_label": hit.get("unit_label", ""),
                            "query_id": hit.get("query_id", ""),
                            "call_fraction": hit.get("call_fraction", 0.0),
                            "diagnostic_matched": matched,
                            "diagnostic_total": total,
                            "diagnostic_missing": missing_count(total, matched),
                            "diagnostic_fraction": fraction(matched, total),
                            "first_pass_distinct": hit.get("first_pass_distinct", 0),
                            "first_pass_diagnostic_total": hit.get("first_pass_diagnostic_total", 0),
                            "member_genes": stats.get("member_genes", ""),
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
                "min_gene_threshold": gene_thr,
                "min_report_unit_threshold": unit_thr,
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
    write_csv(args.out_dir / "truth_kmer_evidence.csv", list(truth_evidence_rows[0].keys()) if truth_evidence_rows else [], truth_evidence_rows)

    missed_counts = collections.Counter(
        (
            row["mode"],
            row["k"],
            row["min_gene_threshold"],
            row["min_report_unit_threshold"],
            row["report_unit"],
        )
        for row in missed_unit_rows
    )
    missed_priority_rows = sorted(
        missed_unit_rows,
        key=lambda row: (
            missed_counts[
                (
                    row["mode"],
                    row["k"],
                    row["min_gene_threshold"],
                    row["min_report_unit_threshold"],
                    row["report_unit"],
                )
            ],
            numeric(row.get("best_diagnostic_fraction", row.get("diagnostic_fraction", 0))),
            numeric(row.get("best_diagnostic_matched", row.get("diagnostic_matched", 0))),
            numeric(row.get("identity", 0)),
            numeric(row.get("coverage", 0)),
        ),
        reverse=True,
    )
    write_csv(args.out_dir / "missed_truth_priority.csv", list(missed_priority_rows[0].keys()) if missed_priority_rows else [], missed_priority_rows)

    missed_summary_counter: dict[tuple[str, str, str, str, str], collections.Counter[str]] = collections.defaultdict(collections.Counter)
    for row in missed_unit_rows:
        key = (
            str(row.get("min_gene_threshold", "")),
            str(row.get("min_report_unit_threshold", "")),
            row.get("report_unit", ""),
            row.get("method", ""),
            row.get("recall_failure_category", ""),
        )
        missed_summary_counter[key]["count"] += 1
    missed_summary_rows = [
        {
            "min_gene_threshold": gene_threshold,
            "min_report_unit_threshold": unit_threshold,
            "report_unit": report_unit,
            "method": method,
            "recall_failure_category": category,
            "misses": counts["count"],
        }
        for (gene_threshold, unit_threshold, report_unit, method, category), counts in sorted(
            missed_summary_counter.items()
        )
    ]
    write_csv(args.out_dir / "missed_truth_summary.csv", list(missed_summary_rows[0].keys()) if missed_summary_rows else [], missed_summary_rows)

    per_unit_rows = [
        {
            "mode": mode,
            "k": k,
            "min_gene_threshold": gene_threshold,
            "min_report_unit_threshold": unit_threshold,
            "report_unit": unit,
            "tp": counts["tp"],
            "fp": counts["fp"],
            "fn": counts["fn"],
        }
        for (mode, k, gene_threshold, unit_threshold, unit), counts in sorted(per_unit.items())
    ]
    write_csv(args.out_dir / "report_unit_failure_summary.csv", list(per_unit_rows[0].keys()) if per_unit_rows else [], per_unit_rows)


if __name__ == "__main__":
    main()
