from __future__ import annotations

import argparse
from pathlib import Path

try:
    from .common import ensure_dir, read_csv
except ImportError:
    from common import ensure_dir, read_csv


def main() -> None:
    parser = argparse.ArgumentParser(description="Create Markdown summary from AMR benchmark outputs")
    parser.add_argument("--selected-manifest", type=Path, required=True)
    parser.add_argument("--aggregate-metrics", type=Path, required=True)
    parser.add_argument("--species-metrics", type=Path)
    parser.add_argument("--out-md", type=Path, required=True)
    args = parser.parse_args()

    selected = read_csv(args.selected_manifest)
    aggregate = read_csv(args.aggregate_metrics)
    species = read_csv(args.species_metrics) if args.species_metrics and args.species_metrics.exists() else []
    if not aggregate:
        raise SystemExit(f"No aggregate rows found in {args.aggregate_metrics}")
    best_exact = max(aggregate, key=lambda row: float(row["exact_f1"]))
    best_report_unit = max(aggregate, key=lambda row: float(row["report_unit_f1"]))

    def describe(row: dict[str, str]) -> str:
        if row.get("threshold_mode") == "fraction":
            return (
                f"mode={row.get('mode', 'direct')}, k={row['k']}, min_gene_fraction={row['min_gene_threshold']}, "
                f"min_report_unit_fraction={row.get('min_report_unit_threshold', row.get('min_gene_group_threshold', ''))}"
            )
        return (
            f"mode={row.get('mode', 'direct')}, k={row['k']}, min_gene_hits={row['min_gene_threshold']}, "
            f"min_report_unit_hits={row.get('min_report_unit_threshold', row.get('min_gene_group_threshold', ''))}"
        )

    species_counts = {}
    for row in selected:
        species_counts[row["species"]] = species_counts.get(row["species"], 0) + 1
    top_species = sorted(species_counts.items(), key=lambda item: (-item[1], item[0]))[:15]

    lines = [
        "# AMR Benchmark Summary",
        "",
        f"- Selected assemblies: {len(selected)}",
        f"- Best exact-F1 config: {describe(best_exact)}, exact_f1={float(best_exact['exact_f1']):.4f}, exact_sensitivity={float(best_exact['exact_sensitivity']):.4f}, exact_specificity={float(best_exact['exact_specificity']):.4f}",
        f"- Best report-unit-F1 config: {describe(best_report_unit)}, report_unit_f1={float(best_report_unit['report_unit_f1']):.4f}, report_unit_sensitivity={float(best_report_unit['report_unit_sensitivity']):.4f}, report_unit_specificity={float(best_report_unit['report_unit_specificity']):.4f}",
        "",
        "## Top Species in Selected Cohort",
        "",
    ]
    for species_name, count in top_species:
        lines.append(f"- {species_name}: {count}")

    lines.extend(
        [
            "",
            "## Parameter Summary",
            "",
            "| mode | k | gene_threshold | report_unit_threshold | exact_f1 | exact_sensitivity | exact_specificity | report_unit_f1 | report_unit_sensitivity | report_unit_specificity | assemblies_compared |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in sorted(
        aggregate,
        key=lambda item: (
            item.get("mode", "direct"),
            item["threshold_mode"],
            int(item["k"]),
            float(item["min_gene_threshold"]),
            float(item.get("min_report_unit_threshold", item.get("min_gene_group_threshold", "0"))),
        ),
    ):
        lines.append(
            f"| {row.get('mode', 'direct')}:{row['threshold_mode']} | {row['k']} | {row['min_gene_threshold']} | {row.get('min_report_unit_threshold', row.get('min_gene_group_threshold', ''))} | {float(row['exact_f1']):.4f} | {float(row['exact_sensitivity']):.4f} | {float(row['exact_specificity']):.4f} | {float(row['report_unit_f1']):.4f} | {float(row['report_unit_sensitivity']):.4f} | {float(row['report_unit_specificity']):.4f} | {row['assemblies_compared']} |"
        )

    ensure_dir(args.out_md.parent)
    args.out_md.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
