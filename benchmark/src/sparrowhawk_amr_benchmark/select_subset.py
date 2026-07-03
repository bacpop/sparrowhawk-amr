from __future__ import annotations

import argparse
import sys
import collections
import hashlib
import math
import random
from pathlib import Path

try:
    from .common import (
        AssemblyRecord,
        MANIFEST_COLUMNS,
        assembly_manifest_row,
        json_dump,
        load_assembly_records,
        quantile_bin,
        write_csv,
    )
except ImportError:
    from common import (  # type: ignore
        AssemblyRecord,
        MANIFEST_COLUMNS,
        assembly_manifest_row,
        json_dump,
        load_assembly_records,
        quantile_bin,
        write_csv,
    )


DEFAULT_ESKAPEE_TARGETS = [
    "Enterococcus faecium",
    "Streptococcus pneumoniae",
    "Klebsiella pneumoniae",
    "Acinetobacter baumannii",
    "Pseudomonas aeruginosa",
    "Enterobacter spp.",
    "Escherichia coli",
]

SPECIES_ALIASES = {
    "Acinetobacter baumanii": "Acinetobacter baumannii",
    "Pseudomonas aeroginosa": "Pseudomonas aeruginosa",
}


def allocate_quotas(
    groups: dict[str, list[AssemblyRecord]],
    target: int,
    min_per_species: int,
    max_per_species: int,
) -> dict[str, int]:
    species = sorted(groups)
    quotas = {sp: min(len(groups[sp]), min_per_species) for sp in species}
    assigned = sum(quotas.values())
    effective_target = max(target, assigned)

    remaining = effective_target - assigned
    if remaining <= 0:
        return quotas

    weights = {sp: math.sqrt(len(groups[sp])) for sp in species}
    weight_total = sum(weights.values()) or 1.0
    fractional = []
    for sp in species:
        headroom = min(len(groups[sp]), max_per_species) - quotas[sp]
        if headroom <= 0:
            continue
        raw = remaining * (weights[sp] / weight_total)
        extra = min(headroom, int(math.floor(raw)))
        quotas[sp] += extra
        fractional.append((raw - extra, sp))

    assigned = sum(quotas.values())
    leftover = effective_target - assigned
    for _, sp in sorted(fractional, reverse=True):
        if leftover <= 0:
            break
        headroom = min(len(groups[sp]), max_per_species) - quotas[sp]
        if headroom <= 0:
            continue
        quotas[sp] += 1
        leftover -= 1

    while leftover > 0:
        progressed = False
        for sp in sorted(species, key=lambda name: (len(groups[name]), name), reverse=True):
            headroom = min(len(groups[sp]), max_per_species) - quotas[sp]
            if headroom <= 0:
                continue
            quotas[sp] += 1
            leftover -= 1
            progressed = True
            if leftover <= 0:
                break
        if not progressed:
            break
    return quotas


def choose_within_species(records: list[AssemblyRecord], quota: int, seed: int) -> list[AssemblyRecord]:
    rng = random.Random(seed)
    bins = quantile_bin(records, bins=4)
    for bucket in bins:
        rng.shuffle(bucket)
        bucket.sort(key=lambda record: (-record.richness_score, -record.n_genes, record.assembly_id))
    chosen: list[AssemblyRecord] = []
    idx = 0
    while len(chosen) < quota:
        progressed = False
        for bucket in bins:
            if idx < len(bucket) and len(chosen) < quota:
                chosen.append(bucket[idx])
                progressed = True
        if not progressed:
            break
        idx += 1
    if len(chosen) < quota:
        seen = {record.assembly_id for record in chosen}
        for record in sorted(records, key=lambda r: (-r.richness_score, -r.n_genes, r.assembly_id)):
            if record.assembly_id in seen:
                continue
            chosen.append(record)
            if len(chosen) >= quota:
                break
    return chosen[:quota]


def stable_species_seed(base_seed: int, species: str) -> int:
    digest = hashlib.sha256(species.encode("utf-8")).hexdigest()
    return base_seed + int(digest[:8], 16)


def split_csv(raw: str) -> list[str]:
    return [part.strip() for part in raw.split(",") if part.strip()]


def canonical_species_name(value: str) -> str:
    return SPECIES_ALIASES.get(value, value)


def matches_species_target(record: AssemblyRecord, target: str) -> bool:
    target = canonical_species_name(target)
    species = canonical_species_name(record.species)
    if target == "Enterobacter spp.":
        return record.genus == "Enterobacter" or species.startswith("Enterobacter ")
    return species == target


def classes_in_records(records: list[AssemblyRecord]) -> list[str]:
    return sorted({class_name for record in records for class_name in record.classes if class_name})


def ranked_records(records: list[AssemblyRecord], seed: int) -> list[AssemblyRecord]:
    rng = random.Random(seed)
    shuffled = records[:]
    rng.shuffle(shuffled)
    return sorted(shuffled, key=lambda record: (-record.richness_score, -record.n_genes, record.assembly_id))


def add_records_for_floor(
    selected: list[AssemblyRecord],
    candidates: list[AssemblyRecord],
    required_count: int,
    seed: int,
) -> None:
    selected_ids = {record.assembly_id for record in selected}
    current = sum(1 for record in selected if record.assembly_id in {candidate.assembly_id for candidate in candidates})
    if current >= required_count:
        return
    for record in ranked_records(candidates, seed):
        if record.assembly_id in selected_ids:
            continue
        selected.append(record)
        selected_ids.add(record.assembly_id)
        current += 1
        if current >= required_count:
            break


def floor_availability(
    records: list[AssemblyRecord],
    species_targets: list[str],
    species_floor: int,
    class_targets: list[str],
    class_floor: int,
) -> tuple[dict[str, int], dict[str, int], list[dict[str, object]]]:
    species_effective = {}
    class_effective = {}
    shortfalls: list[dict[str, object]] = []

    for target in species_targets:
        available = sum(1 for record in records if matches_species_target(record, target))
        effective = min(species_floor, available)
        species_effective[target] = effective
        if available < species_floor:
            shortfalls.append(
                {
                    "kind": "species",
                    "name": target,
                    "required": species_floor,
                    "available": available,
                    "effective_required": effective,
                    "selected": 0,
                }
            )

    for class_name in class_targets:
        available = sum(1 for record in records if class_name in record.classes)
        effective = min(class_floor, available)
        class_effective[class_name] = effective
        if available < class_floor:
            shortfalls.append(
                {
                    "kind": "class",
                    "name": class_name,
                    "required": class_floor,
                    "available": available,
                    "effective_required": effective,
                    "selected": 0,
                }
            )

    return species_effective, class_effective, shortfalls


def add_records_to_target(
    selected: list[AssemblyRecord],
    records: list[AssemblyRecord],
    target_size: int,
    seed: int,
) -> None:
    selected_ids = {record.assembly_id for record in selected}
    if len(selected_ids) >= target_size:
        return
    for record in ranked_records(records, seed):
        if record.assembly_id in selected_ids:
            continue
        selected.append(record)
        selected_ids.add(record.assembly_id)
        if len(selected_ids) >= target_size:
            break


def select_records(args: argparse.Namespace, records: list[AssemblyRecord]) -> tuple[list[AssemblyRecord], dict[str, object]]:
    if args.full_set:
        selected = sorted(records, key=lambda record: (record.species, -record.richness_score, record.assembly_id))
        return selected, {"selection_mode": "full_set", "target_size": len(selected)}

    species_targets = split_csv(args.eskapee_species) if args.eskapee_species else []
    class_targets = split_csv(args.antibiotic_classes) if args.antibiotic_classes else classes_in_records(records)
    species_effective_floors, class_effective_floors, shortfalls = floor_availability(
        records,
        species_targets,
        args.eskapee_floor,
        class_targets,
        args.class_floor,
    )

    groups: dict[str, list[AssemblyRecord]] = collections.defaultdict(list)
    for record in records:
        groups[record.species].append(record)

    quotas = allocate_quotas(groups, args.target_size, args.min_per_species, args.max_per_species)
    selected: list[AssemblyRecord] = []
    for species, quota in sorted(quotas.items()):
        if quota <= 0:
            continue
        selected.extend(choose_within_species(groups[species], quota, stable_species_seed(args.seed, species)))

    for target, floor in species_effective_floors.items():
        candidates = [record for record in records if matches_species_target(record, target)]
        add_records_for_floor(selected, candidates, floor, stable_species_seed(args.seed, target))

    for class_name, floor in class_effective_floors.items():
        candidates = [record for record in records if class_name in record.classes]
        add_records_for_floor(selected, candidates, floor, stable_species_seed(args.seed, class_name))

    deduped = {record.assembly_id: record for record in selected}
    selected = sorted(deduped.values(), key=lambda record: (record.species, -record.richness_score, record.assembly_id))
    pre_topup_size = len(selected)
    add_records_to_target(selected, records, args.target_size, args.seed)
    deduped = {record.assembly_id: record for record in selected}
    selected = sorted(deduped.values(), key=lambda record: (record.species, -record.richness_score, record.assembly_id))
    topup_count = len(selected) - pre_topup_size
    for item in shortfalls:
        if item["kind"] == "species":
            selected_count = sum(1 for record in selected if matches_species_target(record, str(item["name"])))
        else:
            selected_count = sum(1 for record in selected if str(item["name"]) in record.classes)
        item["selected"] = selected_count
        print(
            "Warning: unsatisfied floor "
            f"{item['kind']}={item['name']} required={item['required']} available={item['available']}; "
            f"using_effective_required={item['effective_required']}. "
            "The remaining target size is filled from other available assemblies.",
            file=sys.stderr,
        )
    if len(selected) < args.target_size:
        print(
            f"Warning: requested target_size={args.target_size}, but only {len(selected)} "
            "unique assemblies are available after feasible floor selection and top-up.",
            file=sys.stderr,
        )
    return selected, {
        "selection_mode": "stratified",
        "target_size": args.target_size,
        "pre_topup_size": pre_topup_size,
        "topup_count": topup_count,
        "target_reached": len(selected) >= args.target_size,
        "effective_eskapee_floors": species_effective_floors,
        "effective_class_floors": class_effective_floors,
        "species_quotas": quotas,
        "eskapee_species": species_targets,
        "eskapee_floor": args.eskapee_floor,
        "antibiotic_classes": class_targets,
        "class_floor": args.class_floor,
        "shortfalls": shortfalls,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Select AMR benchmark cohort")
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--target-size", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--min-per-species", type=int, default=25)
    parser.add_argument("--max-per-species", type=int, default=300)
    parser.add_argument("--full-set", action="store_true", help="Use every assembly in the input CSV")
    parser.add_argument("--eskapee-species", default=",".join(DEFAULT_ESKAPEE_TARGETS))
    parser.add_argument("--eskapee-floor", type=int, default=300)
    parser.add_argument("--antibiotic-classes", default="", help="Comma-separated class names; default is all classes in the CSV")
    parser.add_argument("--class-floor", type=int, default=300)
    args = parser.parse_args()

    records = load_assembly_records(args.csv)
    selected, summary = select_records(args, records)
    selected_ids = {record.assembly_id for record in selected}
    reserve = [
        record
        for record in sorted(records, key=lambda r: (r.species, -r.richness_score, r.assembly_id))
        if record.assembly_id not in selected_ids
    ][: max(200, args.target_size // 10)]

    write_csv(
        args.out_csv,
        MANIFEST_COLUMNS,
        (assembly_manifest_row(record, idx + 1) for idx, record in enumerate(selected)),
    )

    json_dump(
        args.out_json,
        {
            **summary,
            "input_assemblies": len(records),
            "selected_size": len(selected),
            "reserve_size": len(reserve),
            "top_species_selected": collections.Counter(record.species for record in selected).most_common(20),
            "class_counts_selected": collections.Counter(class_name for record in selected for class_name in record.classes).most_common(),
            "reserve_assemblies": [
                {
                    "assembly_id": record.assembly_id,
                    "species": record.species,
                    "richness_score": record.richness_score,
                }
                for record in reserve
            ],
        },
    )


if __name__ == "__main__":
    main()
