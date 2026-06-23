from __future__ import annotations

import argparse
import collections
from pathlib import Path

try:
    from .common import AssemblyRecord, json_dump, load_assembly_records, write_csv
except ImportError:
    from common import AssemblyRecord, json_dump, load_assembly_records, write_csv


def build_summary(records: list[AssemblyRecord]) -> dict:
    species_counter = collections.Counter(record.species for record in records)
    genus_counter = collections.Counter(record.genus for record in records)
    rows_by_species = collections.Counter()
    rows_by_class = collections.Counter()
    for record in records:
        rows_by_species[record.species] += record.row_count
        for class_name in record.classes:
            rows_by_class[class_name] += 1
    return {
        "unique_assemblies": len(records),
        "top_species_by_assemblies": species_counter.most_common(20),
        "top_genera_by_assemblies": genus_counter.most_common(20),
        "top_species_by_rows": rows_by_species.most_common(20),
        "top_classes_by_assemblies": rows_by_class.most_common(20),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Profile AMR CSV at assembly level")
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    args = parser.parse_args()

    records = load_assembly_records(args.csv)
    summary = build_summary(records)

    json_dump(args.out_json, summary)
    write_csv(
        args.out_csv,
        [
            "assembly_id",
            "genus",
            "species",
            "organism",
            "taxon_id",
            "row_count",
            "n_classes",
            "n_genes",
            "richness_score",
            "classes",
            "genes",
            "antibiotics",
            "databases",
        ],
        (
            {
                "assembly_id": record.assembly_id,
                "genus": record.genus,
                "species": record.species,
                "organism": record.organism,
                "taxon_id": record.taxon_id,
                "row_count": record.row_count,
                "n_classes": record.n_classes,
                "n_genes": record.n_genes,
                "richness_score": f"{record.richness_score:.3f}",
                "classes": ";".join(record.classes),
                "genes": ";".join(record.genes),
                "antibiotics": ";".join(record.antibiotics),
                "databases": ";".join(record.databases),
            }
            for record in records
        ),
    )


if __name__ == "__main__":
    main()
