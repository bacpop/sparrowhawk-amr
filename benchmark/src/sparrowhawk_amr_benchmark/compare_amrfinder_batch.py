from __future__ import annotations

import argparse
import collections
import csv
import json
import re
from pathlib import Path
from typing import Any

try:
    from .common import f1, read_csv, safe_ratio, write_csv
except ImportError:
    from common import f1, read_csv, safe_ratio, write_csv


OLD_STATUS_RE = re.compile(r"mg_(\d+)__mf_(\d+)_status$")
NEW_STATUS_RE = re.compile(r"gf_(\d+p\d{3})__ff_(\d+p\d{3})_status$")
METRICS = ("exact", "report_unit")


def parse_fraction_label(raw: str) -> float:
    whole, frac = raw.split("p", 1)
    return int(whole) + (int(frac) / 1000.0)


def parse_status_name(path: Path, detector_root: Path) -> dict[str, object]:
    try:
        relative_parts = path.relative_to(detector_root).parts
    except ValueError as exc:
        raise ValueError(
            f"Detector status CSV must be under --detector-root {detector_root}; got {path}. "
            "Pass the detector *_status.csv files, not the comparison output status.csv."
        ) from exc
    mode = "direct"
    k_part = path.parent.name
    if len(relative_parts) >= 3 and not relative_parts[0].startswith("k_"):
        mode = relative_parts[0]
        k_part = relative_parts[1]
    k_match = re.fullmatch(r"k_(\d+)", k_part)
    if not k_match:
        raise ValueError(f"Unexpected detector status filename: {path}")
    stem = path.stem
    old_match = OLD_STATUS_RE.fullmatch(stem)
    if old_match:
        return {
            "mode": mode,
            "k": int(k_match.group(1)),
            "threshold_mode": "absolute",
            "min_gene_threshold": old_match.group(1),
            "min_report_unit_threshold": old_match.group(2),
            "min_gene_fraction": "",
            "min_report_unit_fraction": "",
        }
    new_match = NEW_STATUS_RE.fullmatch(stem)
    if new_match:
        min_gene_fraction = parse_fraction_label(new_match.group(1))
        min_report_unit_fraction = parse_fraction_label(new_match.group(2))
        return {
            "mode": mode,
            "k": int(k_match.group(1)),
            "threshold_mode": "fraction",
            "min_gene_threshold": f"{min_gene_fraction:.3f}",
            "min_report_unit_threshold": f"{min_report_unit_fraction:.3f}",
            "min_gene_fraction": f"{min_gene_fraction:.3f}",
            "min_report_unit_fraction": f"{min_report_unit_fraction:.3f}",
        }
    raise ValueError(f"Unexpected detector status filename: {path}")


def load_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def resolve_existing_path(raw: str, base_dir: Path) -> Path:
    path = Path(raw)
    if path.is_absolute():
        return path
    candidate = base_dir / path
    if candidate.exists():
        return candidate
    return path


def reportable_rows(path: Path, included_types: set[str]) -> list[dict[str, str]]:
    return [
        row
        for row in load_tsv(path)
        if row.get("Type", "").strip().upper() in included_types
    ]


def load_hierarchy(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    rows = load_tsv(path)
    return {row["node_id"]: row for row in rows if row.get("node_id")}


def gene_group_key(node: str, hierarchy: dict[str, dict[str, str]]) -> str:
    node = node.strip()
    if not node:
        return ""
    meta = hierarchy.get(node)
    if meta and meta.get("allele") == "1" and meta.get("parent_node_id"):
        return meta["parent_node_id"]
    return node


def load_report_map_rows(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def report_map_lookup(rows: list[dict[str, str]]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for row in rows:
        accession = row.get("protein_accession", "")
        symbol = row.get("element_symbol", "")
        node = row.get("hierarchy_node", "")
        unit = row.get("report_unit_key", "")
        for key in (accession, symbol, node):
            if key and unit:
                mapping.setdefault(key, unit)
    return mapping


def metric_universes(
    report_map_rows: list[dict[str, str]],
    hierarchy: dict[str, dict[str, str]],
    included_types: set[str],
) -> dict[str, set[str]]:
    report_map_rows = [
        row
        for row in report_map_rows
        if row.get("type", "").strip().upper() in included_types
    ]
    exact = {
        row.get("element_symbol", "")
        for row in report_map_rows
        if row.get("element_symbol", "")
    }
    report_unit = {
        row.get("report_unit_key", "")
        for row in report_map_rows
        if row.get("report_unit_key", "")
    }
    return {
        "exact": {value for value in exact if value},
        "report_unit": {value for value in report_unit if value},
    }


def report_map_path(report_map_root: Path, mode: str, k: int) -> Path:
    alphabet = "protein" if mode == "protein_cds" else "dna"
    return report_map_root / f"{alphabet}_k{k}.tsv"


def normalize_amrfinder(
    rows: list[dict[str, str]],
    report_map: dict[str, str],
    hierarchy: dict[str, dict[str, str]],
) -> dict[str, list[str]]:
    exact = sorted({row["Element symbol"] for row in rows if row.get("Element symbol")})
    gene_groups = sorted(
        {
            gene_group_key(row.get("Hierarchy node") or row.get("Element symbol", ""), hierarchy)
            for row in rows
            if row.get("Hierarchy node") or row.get("Element symbol")
        }
    )
    report_units = set()
    for row in rows:
        accession = row.get("Closest reference accession", "")
        symbol = row.get("Element symbol", "")
        node = row.get("Hierarchy node", "")
        unit = report_map.get(accession) or report_map.get(symbol) or report_map.get(node)
        if unit:
            report_units.add(unit)
        elif node:
            report_units.add(f"hierarchy_node:{gene_group_key(node, hierarchy)}")
        elif symbol:
            report_units.add(f"exact_gene:{symbol}")
    return {
        "exact": exact,
        "gene_group": sorted(gene_groups),
        "report_unit": sorted(report_units),
    }


def detector_hit_type(hit: dict[str, Any]) -> str:
    return str(hit.get("type_name") or hit.get("type") or "").strip().upper()


def filter_detector_hits_by_type(payload: dict[str, Any], included_types: set[str]) -> dict[str, Any]:
    return {
        **payload,
        "hits": [
            hit
            for hit in payload.get("hits", [])
            if detector_hit_type(hit) in included_types
        ],
    }


def normalize_detector(
    payload: dict[str, Any],
    hierarchy: dict[str, dict[str, str]],
) -> dict[str, list[str]]:
    exact = set()
    gene_groups = set()
    report_units = set()
    for hit in payload.get("hits", []):
        if hit.get("element_symbol"):
            exact.add(hit["element_symbol"])
        elif hit.get("gene_id") and (hit.get("unit_type") == "exact_gene" or hit.get("call_type") == "gene"):
            exact.add(hit["gene_id"].split("|", 1)[0])

        node = hit.get("hierarchy_node") or hit.get("gene_group") or hit.get("element_symbol") or ""
        key = gene_group_key(node, hierarchy)
        if key:
            gene_groups.add(key)

        if hit.get("unit_id"):
            if hit.get("unit_type"):
                report_units.add(f"{hit['unit_type']}:{hit['unit_id']}")
            elif hit.get("call_type") == "gene":
                report_units.add(f"exact_gene:{hit['unit_id']}")
            elif hit.get("call_type") in {"gene_group", "family"}:
                report_units.add(f"hierarchy_node:{hit['unit_id']}")
    return {
        "exact": sorted(exact),
        "gene_group": sorted(gene_groups),
        "report_unit": sorted(report_units),
    }


def counts(detector: set[str], baseline: set[str], universe: set[str]) -> dict[str, int]:
    tp = len(detector & baseline)
    fp = len(detector - baseline)
    fn = len(baseline - detector)
    tn = len(universe - (detector | baseline))
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "pred": len(detector), "truth": len(baseline)}


def add_metric_fields(prefix: str, values: collections.Counter[str]) -> dict[str, object]:
    tp = values[f"{prefix}_tp"]
    fp = values[f"{prefix}_fp"]
    fn = values[f"{prefix}_fn"]
    tn = values[f"{prefix}_tn"]
    pred = values[f"{prefix}_pred"]
    truth = values[f"{prefix}_truth"]
    return {
        f"{prefix}_precision": safe_ratio(tp, pred),
        f"{prefix}_recall": safe_ratio(tp, truth),
        f"{prefix}_sensitivity": safe_ratio(tp, tp + fn),
        f"{prefix}_specificity": safe_ratio(tn, tn + fp),
        f"{prefix}_f1": f1(tp, pred, truth),
        f"{prefix}_tp": tp,
        f"{prefix}_fp": fp,
        f"{prefix}_fn": fn,
        f"{prefix}_tn": tn,
        f"{prefix}_pred": pred,
        f"{prefix}_truth": truth,
    }


def add_counts(counter: collections.Counter[str], row_counts: dict[str, dict[str, int]]) -> None:
    counter["assemblies"] += 1
    for metric, metric_counts in row_counts.items():
        for key, value in metric_counts.items():
            counter[f"{metric}_{key}"] += value


def metric_row(base: dict[str, object], counter: collections.Counter[str]) -> dict[str, object]:
    row = {
        **base,
        "assemblies_compared": counter["assemblies"],
    }
    for metric in METRICS:
        row.update(add_metric_fields(metric, counter))
    return row


def empty_universes() -> dict[str, set[str]]:
    return {metric: set() for metric in METRICS}


def label_value(value: object) -> str:
    text = str(value or "").strip()
    return text if text else "Unclassified"


def labels_from_amrfinder_rows(rows: list[dict[str, str]], column: str, include_blank: bool = False) -> set[str]:
    labels = {label_value(row.get(column, "")) for row in rows}
    return labels if include_blank else {label for label in labels if label != "Unclassified"}


def labels_from_detector_hits(payload: dict[str, Any], field: str, include_blank: bool = False) -> set[str]:
    labels = {label_value(hit.get(field)) for hit in payload.get("hits", [])}
    return labels if include_blank else {label for label in labels if label != "Unclassified"}


def normalize_detector_hits(
    hits: list[dict[str, Any]],
    hierarchy: dict[str, dict[str, str]],
) -> dict[str, list[str]]:
    exact = set()
    gene_groups = set()
    report_units = set()
    for hit in hits:
        if hit.get("element_symbol"):
            exact.add(hit["element_symbol"])
        elif hit.get("gene_id") and (hit.get("unit_type") == "exact_gene" or hit.get("call_type") == "gene"):
            exact.add(hit["gene_id"].split("|", 1)[0])

        node = hit.get("hierarchy_node") or hit.get("gene_group") or hit.get("element_symbol") or ""
        key = gene_group_key(node, hierarchy)
        if key:
            gene_groups.add(key)

        if hit.get("unit_id"):
            if hit.get("unit_type"):
                report_units.add(f"{hit['unit_type']}:{hit['unit_id']}")
            elif hit.get("call_type") == "gene":
                report_units.add(f"exact_gene:{hit['unit_id']}")
            elif hit.get("call_type") in {"gene_group", "family"}:
                report_units.add(f"hierarchy_node:{hit['unit_id']}")
    return {
        "exact": sorted(exact),
        "gene_group": sorted(gene_groups),
        "report_unit": sorted(report_units),
    }


def filter_native_rows(rows: list[dict[str, str]], column: str, label: str) -> list[dict[str, str]]:
    return [row for row in rows if label_value(row.get(column, "")) == label]


def filter_detector_hits(payload: dict[str, Any], field: str, label: str) -> list[dict[str, Any]]:
    return [hit for hit in payload.get("hits", []) if label_value(hit.get(field)) == label]


def update_universe(universe: dict[str, set[str]], detector_norm: dict[str, list[str]], baseline_norm: dict[str, list[str]]) -> None:
    for metric in METRICS:
        universe[metric].update(detector_norm[metric])
        universe[metric].update(baseline_norm[metric])


def metric_counts(detector_norm: dict[str, list[str]], baseline_norm: dict[str, list[str]], universes: dict[str, set[str]]) -> dict[str, dict[str, int]]:
    return {
        metric: counts(
            set(detector_norm[metric]),
            set(baseline_norm[metric]),
            universes[metric],
        )
        for metric in METRICS
    }


def grouped_rows(
    grouped: dict[str, collections.Counter[str]],
    params: dict[str, object],
    group_columns: tuple[str, ...],
) -> list[dict[str, object]]:
    rows = []
    for key, counter in sorted(grouped.items()):
        parts = key.split("||")
        base = dict(zip(group_columns, parts))
        base.update(params)
        rows.append(metric_row(base, counter))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare detector outputs to native AMRFinderPlus TSVs")
    parser.add_argument("--amrfinder-status", type=Path, required=True)
    parser.add_argument("--detector-root", type=Path, required=True)
    parser.add_argument("--report-map-root", type=Path, required=True)
    parser.add_argument("--hierarchy", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--status-csv", type=Path, action="append", required=True)
    parser.add_argument("--include-types", default="AMR,STRESS,VIRULENCE")
    args = parser.parse_args()

    included_types = {value.strip().upper() for value in args.include_types.split(",") if value.strip()}
    hierarchy = load_hierarchy(args.hierarchy)
    amrfinder_status = {}
    for row in read_csv(args.amrfinder_status):
        if row.get("returncode") != "0":
            continue
        if row.get("tsv_path"):
            row = {**row, "tsv_path": str(resolve_existing_path(row["tsv_path"], args.amrfinder_status.parent))}
        amrfinder_status[row["assembly_id"]] = row
    status_paths = args.status_csv
    params_by_status = {
        status: parse_status_name(status, args.detector_root)
        for status in status_paths
    }

    report_map_by_key: dict[tuple[str, int], dict[str, str]] = {}
    universes_by_key: dict[tuple[str, int], dict[str, set[str]]] = {}
    for params in params_by_status.values():
        key = (str(params["mode"]), int(params["k"]))
        if key in report_map_by_key:
            continue
        rows = load_report_map_rows(report_map_path(args.report_map_root, key[0], key[1]))
        report_map_by_key[key] = report_map_lookup(rows)
        universes_by_key[key] = metric_universes(rows, hierarchy, included_types)

    aggregate_rows = []
    species_rows = []
    class_rows = []
    subclass_rows = []
    species_class_rows = []
    type_rows = []

    for status_csv in status_paths:
        params = params_by_status[status_csv]
        key = (str(params["mode"]), int(params["k"]))
        report_map = report_map_by_key[key]
        universes = universes_by_key[key]
        micro: collections.Counter[str] = collections.Counter()
        species_micro: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)

        class_items = []
        subclass_items = []
        species_class_items = []
        type_items = []
        class_universes: dict[str, dict[str, set[str]]] = collections.defaultdict(empty_universes)
        subclass_universes: dict[str, dict[str, set[str]]] = collections.defaultdict(empty_universes)
        species_class_universes: dict[str, dict[str, set[str]]] = collections.defaultdict(empty_universes)
        type_universes: dict[str, dict[str, set[str]]] = collections.defaultdict(empty_universes)
        per_assembly = []
        missing_native_class = 0
        missing_detector_class = 0

        for row in read_csv(status_csv):
            assembly_id = row["assembly_id"]
            baseline_row = amrfinder_status.get(assembly_id)
            if not baseline_row or row.get("detector_status") not in {"ok", "cached"}:
                continue

            native_rows = reportable_rows(Path(baseline_row["tsv_path"]), included_types)
            detector_json_path = resolve_existing_path(row["detector_json"], status_csv.parent)
            detector_payload = filter_detector_hits_by_type(load_json(detector_json_path), included_types)
            missing_native_class += sum(1 for native_row in native_rows if not str(native_row.get("Class", "")).strip())
            missing_detector_class += sum(1 for hit in detector_payload.get("hits", []) if not str(hit.get("class_name", "") or "").strip())
            baseline_norm = normalize_amrfinder(native_rows, report_map, hierarchy)
            detector_norm = normalize_detector(detector_payload, hierarchy)
            row_counts = metric_counts(detector_norm, baseline_norm, universes)
            add_counts(micro, row_counts)

            species = row.get("species", "") or "Unknown species"
            add_counts(species_micro[species], row_counts)

            for type_label in labels_from_amrfinder_rows(native_rows, "Type") | labels_from_detector_hits(detector_payload, "type_name"):
                type_baseline = normalize_amrfinder(filter_native_rows(native_rows, "Type", type_label), report_map, hierarchy)
                type_detector = normalize_detector_hits(filter_detector_hits(detector_payload, "type_name", type_label), hierarchy)
                type_items.append((type_label, type_detector, type_baseline))
                update_universe(type_universes[type_label], type_detector, type_baseline)

            for class_label in labels_from_amrfinder_rows(native_rows, "Class") | labels_from_detector_hits(detector_payload, "class_name"):
                class_baseline = normalize_amrfinder(filter_native_rows(native_rows, "Class", class_label), report_map, hierarchy)
                class_detector = normalize_detector_hits(filter_detector_hits(detector_payload, "class_name", class_label), hierarchy)
                class_items.append((class_label, class_detector, class_baseline))
                update_universe(class_universes[class_label], class_detector, class_baseline)

                species_class_label = f"{species}||{class_label}"
                species_class_items.append((species_class_label, class_detector, class_baseline))
                update_universe(species_class_universes[species_class_label], class_detector, class_baseline)

            for subclass_label in labels_from_amrfinder_rows(native_rows, "Subclass") | labels_from_detector_hits(detector_payload, "subclass"):
                subclass_baseline = normalize_amrfinder(filter_native_rows(native_rows, "Subclass", subclass_label), report_map, hierarchy)
                subclass_detector = normalize_detector_hits(filter_detector_hits(detector_payload, "subclass", subclass_label), hierarchy)
                subclass_items.append((subclass_label, subclass_detector, subclass_baseline))
                update_universe(subclass_universes[subclass_label], subclass_detector, subclass_baseline)

            exact_detector = set(detector_norm["exact"])
            exact_baseline = set(baseline_norm["exact"])
            report_detector = set(detector_norm["report_unit"])
            report_baseline = set(baseline_norm["report_unit"])
            per_assembly.append(
                {
                    "assembly_id": assembly_id,
                    "species": row.get("species", ""),
                    "genus": row.get("genus", ""),
                    "classes": row.get("classes", ""),
                    "antibiotics": row.get("antibiotics", ""),
                    **params,
                    "native_amrfinder_tsv": baseline_row["tsv_path"],
                    "detector_json": str(detector_json_path),
                    "baseline_exact_count": len(exact_baseline),
                    "detector_exact_count": len(exact_detector),
                    "exact_tp": row_counts["exact"]["tp"],
                    "exact_fp": row_counts["exact"]["fp"],
                    "exact_fn": row_counts["exact"]["fn"],
                    "exact_tn": row_counts["exact"]["tn"],
                    "exact_sensitivity": safe_ratio(row_counts["exact"]["tp"], row_counts["exact"]["tp"] + row_counts["exact"]["fn"]),
                    "exact_specificity": safe_ratio(row_counts["exact"]["tn"], row_counts["exact"]["tn"] + row_counts["exact"]["fp"]),
                    "report_unit_tp": row_counts["report_unit"]["tp"],
                    "report_unit_fp": row_counts["report_unit"]["fp"],
                    "report_unit_fn": row_counts["report_unit"]["fn"],
                    "report_unit_tn": row_counts["report_unit"]["tn"],
                    "report_unit_sensitivity": safe_ratio(row_counts["report_unit"]["tp"], row_counts["report_unit"]["tp"] + row_counts["report_unit"]["fn"]),
                    "report_unit_specificity": safe_ratio(row_counts["report_unit"]["tn"], row_counts["report_unit"]["tn"] + row_counts["report_unit"]["fp"]),
                    "baseline_only_exact": ";".join(sorted(exact_baseline - exact_detector)),
                    "detector_only_exact": ";".join(sorted(exact_detector - exact_baseline)),
                    "baseline_only_report_unit": ";".join(sorted(report_baseline - report_detector)),
                    "detector_only_report_unit": ";".join(sorted(report_detector - report_baseline)),
                }
            )

        class_micro: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
        for label, detector_norm, baseline_norm in class_items:
            add_counts(class_micro[label], metric_counts(detector_norm, baseline_norm, class_universes[label]))

        subclass_micro: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
        for label, detector_norm, baseline_norm in subclass_items:
            add_counts(subclass_micro[label], metric_counts(detector_norm, baseline_norm, subclass_universes[label]))

        species_class_micro: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
        for label, detector_norm, baseline_norm in species_class_items:
            add_counts(species_class_micro[label], metric_counts(detector_norm, baseline_norm, species_class_universes[label]))

        type_micro: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
        for label, detector_norm, baseline_norm in type_items:
            add_counts(type_micro[label], metric_counts(detector_norm, baseline_norm, type_universes[label]))

        out_name = (
            f"{params['mode']}_k_{params['k']}_"
            f"{params['threshold_mode']}_gene_{params['min_gene_threshold']}_"
            f"report_unit_{params['min_report_unit_threshold']}_assemblies.csv"
        )
        write_csv(args.out_dir / out_name, list(per_assembly[0].keys()) if per_assembly else [], per_assembly)
        if missing_native_class:
            print(f"Warning: {missing_native_class} AMRFinderPlus rows had blank Class values and were omitted from class-level metrics for {status_csv}.")
        if missing_detector_class:
            print(f"Warning: {missing_detector_class} detector hits had blank class_name values and were omitted from class-level metrics for {status_csv}.")
        aggregate_rows.append(metric_row(params, micro))
        species_rows.extend(grouped_rows(species_micro, params, ("species",)))
        class_rows.extend(grouped_rows(class_micro, params, ("class_name",)))
        subclass_rows.extend(grouped_rows(subclass_micro, params, ("subclass",)))
        species_class_rows.extend(grouped_rows(species_class_micro, params, ("species", "class_name")))
        type_rows.extend(grouped_rows(type_micro, params, ("type_name",)))

    write_csv(args.out_dir / "aggregate_metrics.csv", list(aggregate_rows[0].keys()) if aggregate_rows else [], aggregate_rows)
    write_csv(args.out_dir / "species_metrics.csv", list(species_rows[0].keys()) if species_rows else [], species_rows)
    write_csv(args.out_dir / "class_metrics.csv", list(class_rows[0].keys()) if class_rows else [], class_rows)
    write_csv(args.out_dir / "subclass_metrics.csv", list(subclass_rows[0].keys()) if subclass_rows else [], subclass_rows)
    write_csv(args.out_dir / "species_class_metrics.csv", list(species_class_rows[0].keys()) if species_class_rows else [], species_class_rows)
    write_csv(args.out_dir / "type_metrics.csv", list(type_rows[0].keys()) if type_rows else [], type_rows)


if __name__ == "__main__":
    main()
