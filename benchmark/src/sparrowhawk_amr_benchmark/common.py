from __future__ import annotations

import csv
import gzip
import json
import math
import os
import re
import subprocess
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Sequence


CSV_PREFIX = "pheno_geno_merged-"
MANIFEST_COLUMNS = [
    "assembly_id",
    "species",
    "genus",
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
    "selection_rank",
    "selection_group",
]


@dataclass
class AssemblyRecord:
    assembly_id: str
    genus: str
    species: str
    organism: str
    taxon_id: str
    row_count: int
    n_classes: int
    n_genes: int
    richness_score: float
    classes: List[str]
    genes: List[str]
    antibiotics: List[str]
    databases: List[str]


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def parse_multi_source_field(raw: str) -> List[str]:
    if not raw:
        return []
    return sorted({item.strip() for item in raw.split(";") if item.strip()})


def compute_richness_score(
    row_count: int,
    n_classes: int,
    n_genes: int,
    n_antibiotics: int,
    n_databases: int,
) -> float:
    return (
        4.0 * n_classes
        + 2.0 * n_genes
        + 1.5 * n_antibiotics
        + 0.5 * min(row_count, 15)
        + 0.75 * n_databases
    )


def load_assembly_records(csv_path: Path) -> List[AssemblyRecord]:
    assemblies: Dict[str, Dict[str, Any]] = {}
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            assembly_id = row[f"{CSV_PREFIX}assembly_ID"].strip()
            if not assembly_id:
                continue
            record = assemblies.setdefault(
                assembly_id,
                {
                    "assembly_id": assembly_id,
                    "genus": row[f"{CSV_PREFIX}genus"].strip(),
                    "species": row[f"{CSV_PREFIX}species"].strip(),
                    "organism": row[f"{CSV_PREFIX}organism"].strip(),
                    "taxon_id": row[f"{CSV_PREFIX}taxon_id"].strip(),
                    "row_count": 0,
                    "classes": set(),
                    "genes": set(),
                    "antibiotics": set(),
                    "databases": set(),
                },
            )
            record["row_count"] += 1
            if row.get(f"{CSV_PREFIX}class"):
                record["classes"].add(row[f"{CSV_PREFIX}class"].strip())
            if row.get(f"{CSV_PREFIX}gene_symbol"):
                record["genes"].add(row[f"{CSV_PREFIX}gene_symbol"].strip())
            if row.get(f"{CSV_PREFIX}antibiotic_name"):
                record["antibiotics"].add(row[f"{CSV_PREFIX}antibiotic_name"].strip())
            if row.get(f"{CSV_PREFIX}database"):
                record["databases"].update(parse_multi_source_field(row[f"{CSV_PREFIX}database"]))

    out: List[AssemblyRecord] = []
    for record in assemblies.values():
        out.append(
            AssemblyRecord(
                assembly_id=record["assembly_id"],
                genus=record["genus"],
                species=record["species"],
                organism=record["organism"],
                taxon_id=record["taxon_id"],
                row_count=record["row_count"],
                n_classes=len(record["classes"]),
                n_genes=len(record["genes"]),
                richness_score=compute_richness_score(
                    row_count=record["row_count"],
                    n_classes=len(record["classes"]),
                    n_genes=len(record["genes"]),
                    n_antibiotics=len(record["antibiotics"]),
                    n_databases=len(record["databases"]),
                ),
                classes=sorted(record["classes"]),
                genes=sorted(record["genes"]),
                antibiotics=sorted(record["antibiotics"]),
                databases=sorted(record["databases"]),
            )
        )
    return out


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, Any]]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def assembly_manifest_row(record: AssemblyRecord, rank: int, group: str = "primary") -> Dict[str, Any]:
    return {
        "assembly_id": record.assembly_id,
        "species": record.species,
        "genus": record.genus,
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
        "selection_rank": rank,
        "selection_group": group,
    }


def assembly_metadata_by_id(rows: Sequence[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    return {
        row["assembly_id"]: {
            "species": row.get("species", ""),
            "classes": row.get("classes", ""),
            "antibiotics": row.get("antibiotics", ""),
        }
        for row in rows
        if row.get("assembly_id")
    }


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def json_dump(path: Path, payload: Any) -> None:
    ensure_dir(path.parent)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def parse_resfinder_json(path: Path) -> Dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def _collect_gene_hits(obj: Any, out: List[str]) -> None:
    if isinstance(obj, dict):
        if obj.get("gene") is True and isinstance(obj.get("name"), str):
            out.append(obj["name"])
        for value in obj.values():
            _collect_gene_hits(value, out)
    elif isinstance(obj, list):
        for value in obj:
            _collect_gene_hits(value, out)


def gene_family(symbol: str) -> str:
    return symbol.split("_")[0]


def detector_gene_name(symbol: str) -> str:
    return re.sub(r"_\d+(?:_.+)?$", "", symbol)


def normalize_resfinder_hits(payload: Dict[str, Any], report_map: Dict[str, str] | None = None) -> Dict[str, Any]:
    exact: List[str] = []
    _collect_gene_hits(payload, exact)
    exact_set = sorted(set(exact))
    gene_groups = sorted({gene_family(hit) for hit in exact_set})
    report_units = sorted(
        {
            report_map.get(detector_gene_name(hit), f"exact_gene:{detector_gene_name(hit)}")
            if report_map
            else f"exact_gene:{detector_gene_name(hit)}"
            for hit in exact_set
        }
    )
    return {"exact_hits": exact_set, "gene_group_hits": gene_groups, "report_units": report_units}


def normalize_detector_hits(payload: Dict[str, Any]) -> Dict[str, Any]:
    exact = sorted(
        {
            detector_gene_name(hit.get("element_symbol") or hit["gene_id"])
            for hit in payload["hits"]
            if hit.get("gene_id") or hit.get("element_symbol")
        }
    )
    gene_groups = sorted(
        {
            hit.get("gene_family") or hit.get("gene_group")
            for hit in payload["hits"]
            if hit.get("gene_family") or hit.get("gene_group")
        }
    )
    report_units = sorted(
        {
            f"{hit['unit_type']}:{hit['unit_id']}"
            for hit in payload["hits"]
            if hit.get("unit_type") and hit.get("unit_id")
        }
    )
    return {"exact_hits": exact, "gene_group_hits": gene_groups, "report_units": report_units}


def score_overlap(detector: Dict[str, Any], baseline: Dict[str, Any]) -> Dict[str, Any]:
    detector_exact = set(detector["exact_hits"])
    detector_gene_group = set(detector["gene_group_hits"])
    baseline_exact = set(baseline["exact_hits"])
    baseline_gene_group = set(baseline["gene_group_hits"])
    detector_report_units = set(detector.get("report_units", []))
    baseline_report_units = set(baseline.get("report_units", []))

    exact_overlap = sorted(detector_exact & baseline_exact)
    gene_group_overlap = sorted(detector_gene_group & baseline_gene_group)
    report_unit_overlap = sorted(detector_report_units & baseline_report_units)
    detector_only = sorted(detector_exact - baseline_exact)
    baseline_only = sorted(baseline_exact - detector_exact)

    return {
        "detector_hits": sorted(detector_exact),
        "detector_gene_groups": sorted(detector_gene_group),
        "baseline_hits": sorted(baseline_exact),
        "baseline_gene_groups": sorted(baseline_gene_group),
        "detector_report_units": sorted(detector_report_units),
        "baseline_report_units": sorted(baseline_report_units),
        "exact_overlap": exact_overlap,
        "gene_group_overlap": gene_group_overlap,
        "report_unit_overlap": report_unit_overlap,
        "detector_only": detector_only,
        "baseline_only": baseline_only,
        "detector_hit_count": len(detector_exact),
        "detector_gene_group_count": len(detector_gene_group),
        "detector_report_unit_count": len(detector_report_units),
        "baseline_hit_count": len(baseline_exact),
        "baseline_gene_group_count": len(baseline_gene_group),
        "baseline_report_unit_count": len(baseline_report_units),
        "exact_overlap_count": len(exact_overlap),
        "gene_group_overlap_count": len(gene_group_overlap),
        "report_unit_overlap_count": len(report_unit_overlap),
    }


def f1(tp: int, predicted: int, truth: int) -> float:
    if predicted == 0 or truth == 0 or tp == 0:
        return 0.0
    precision = tp / predicted
    recall = tp / truth
    return 2.0 * precision * recall / (precision + recall)


def safe_ratio(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return 0.0
    return numerator / denominator


def shell(command: Sequence[str], cwd: Path | None = None, env: Dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(command),
        cwd=str(cwd) if cwd else None,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


def run_and_time(command: Sequence[str], cwd: Path | None = None, env: Dict[str, str] | None = None) -> Dict[str, Any]:
    start = time.time()
    proc = shell(command, cwd=cwd, env=env)
    elapsed = time.time() - start
    return {
        "returncode": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "elapsed_seconds": round(elapsed, 3),
    }


def fetch_docsum_xml(assembly_accession: str) -> str:
    search = shell(
        [
            "esearch",
            "-db",
            "assembly",
            "-query",
            f"{assembly_accession}[Assembly Accession]",
        ]
    )
    if search.returncode != 0:
        raise RuntimeError(search.stderr.strip() or "esearch failed")
    summary = subprocess.run(
        ["efetch", "-format", "docsum"],
        input=search.stdout,
        text=True,
        capture_output=True,
        check=False,
    )
    if summary.returncode != 0:
        raise RuntimeError(summary.stderr.strip() or "efetch failed")
    return summary.stdout


def resolve_assembly_download_urls(assembly_accession: str) -> list[str]:
    xml_text = fetch_docsum_xml(assembly_accession)
    root = ET.fromstring(xml_text)

    ftp_paths: list[str] = []

    def add_path(value: str | None) -> None:
        text = value.strip() if value else ""
        if text and text not in ftp_paths:
            ftp_paths.append(text)

    for elem in root.iter("FtpPath_RefSeq"):
        add_path(elem.text)
    for elem in root.iter("FtpPath_GenBank"):
        add_path(elem.text)
    for elem in root.iter("FtpPath"):
        add_path(elem.text)

    if not ftp_paths:
        raise RuntimeError(f"No FTP path found for {assembly_accession}")

    urls = []
    for ftp_path in ftp_paths:
        basename = ftp_path.rsplit("/", 1)[-1]
        urls.append(ftp_path.replace("ftp://", "https://") + f"/{basename}_genomic.fna.gz")
    return urls


def resolve_assembly_download_url(assembly_accession: str) -> str:
    return resolve_assembly_download_urls(assembly_accession)[0]


def download_to_path(url: str, destination: Path) -> None:
    ensure_dir(destination.parent)
    result = shell([
        "curl",
        "-L",
        "-f",
        "--retry",
        "3",
        "--retry-delay",
        "5",
        "--retry-connrefused",
        "-o",
        str(destination),
        url,
    ])
    if result.returncode == 0:
        return
    fallback = shell(["wget", "--tries=3", "--waitretry=5", "-O", str(destination), url])
    if fallback.returncode != 0:
        raise RuntimeError(fallback.stderr.strip() or result.stderr.strip() or f"failed to download {url}")


def decompress_gzip(source: Path, destination: Path) -> None:
    ensure_dir(destination.parent)
    with gzip.open(source, "rb") as in_handle, destination.open("wb") as out_handle:
        out_handle.write(in_handle.read())



def quantile_bin(values: Sequence[AssemblyRecord], bins: int = 4) -> List[List[AssemblyRecord]]:
    if not values:
        return []
    ordered = sorted(
        values,
        key=lambda record: (-record.richness_score, -record.n_genes, -record.n_classes, record.assembly_id),
    )
    out: List[List[AssemblyRecord]] = []
    for idx in range(bins):
        start = math.floor(idx * len(ordered) / bins)
        end = math.floor((idx + 1) * len(ordered) / bins)
        out.append(ordered[start:end])
    return [bucket for bucket in out if bucket]
