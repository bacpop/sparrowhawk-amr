from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from sparrowhawk_amr_benchmark.analyze_amrfinder_failures import parse_mode_file
from sparrowhawk_amr_benchmark.compare_amrfinder_batch import filter_detector_hits_by_type
from sparrowhawk_amr_benchmark.run_detector_batch import run_one


class ParseModeFileTests(unittest.TestCase):
    def test_parses_full_threshold_combination(self) -> None:
        path = Path("direct_k_31_fraction_gene_0.100_report_unit_0.100_assemblies.csv")
        self.assertEqual(
            parse_mode_file(path), ("direct", 31, "fraction", "0.100", "0.100")
        )

    def test_parses_protein_cds_mode(self) -> None:
        path = Path("protein_cds_k_5_fraction_gene_0.020_report_unit_0.300_assemblies.csv")
        self.assertEqual(
            parse_mode_file(path), ("protein_cds", 5, "fraction", "0.020", "0.300")
        )

    def test_parses_legacy_gene_group_naming(self) -> None:
        path = Path("cds_k_21_absolute_gene_8_gene_group_12_assemblies.csv")
        self.assertEqual(parse_mode_file(path), ("cds", 21, "absolute", "8", "12"))

    def test_rejects_unexpected_names(self) -> None:
        with self.assertRaises(ValueError):
            parse_mode_file(Path("status.csv"))


class FilterDetectorHitsByTypeTests(unittest.TestCase):
    def test_filters_by_included_types(self) -> None:
        payload = {
            "sample_name": "s1",
            "hits": [
                {"unit_id": "a", "type_name": "AMR"},
                {"unit_id": "b", "type_name": "STRESS"},
            ],
        }
        filtered = filter_detector_hits_by_type(payload, {"AMR"})
        self.assertEqual([hit["unit_id"] for hit in filtered["hits"]], ["a"])

    def test_raises_on_hit_without_type_name(self) -> None:
        payload = {
            "sample_name": "s1",
            "hits": [{"unit_id": "a", "query_id": "contig", "type_name": None}],
        }
        with self.assertRaises(ValueError) as ctx:
            filter_detector_hits_by_type(payload, {"AMR"})
        self.assertIn("without type_name", str(ctx.exception))


class DetectorCacheFingerprintTests(unittest.TestCase):
    def _row(self, fasta_path: Path) -> dict[str, str]:
        return {
            "assembly_id": "GCF_TEST",
            "species": "Testus testus",
            "fetch_status": "cached",
            "local_fasta_path": str(fasta_path),
        }

    def test_cache_hit_requires_matching_index_fingerprint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            result_dir = Path(tmp)
            fasta = result_dir / "GCF_TEST.fna"
            fasta.write_text(">c\nACGT\n")
            (result_dir / "GCF_TEST.json").write_text("{}")
            (result_dir / "GCF_TEST.meta.json").write_text(
                json.dumps({"index_sha256": "sha-current"})
            )
            result = run_one(
                self._row(fasta),
                Path("/nonexistent/detector"),
                Path("/nonexistent/index"),
                "sha-current",
                result_dir,
                0.1,
                0.1,
                "direct",
                result_dir,
            )
            self.assertEqual(result["detector_status"], "cached")

    def test_stale_fingerprint_triggers_rerun(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            result_dir = Path(tmp)
            fasta = result_dir / "GCF_TEST.fna"
            fasta.write_text(">c\nACGT\n")
            (result_dir / "GCF_TEST.json").write_text("{}")
            (result_dir / "GCF_TEST.meta.json").write_text(
                json.dumps({"index_sha256": "sha-old"})
            )
            # The stale cache must not be trusted: run_one goes on to invoke the
            # (deliberately nonexistent) detector binary instead of returning "cached".
            with self.assertRaises(FileNotFoundError):
                run_one(
                    self._row(fasta),
                    Path("/nonexistent/detector"),
                    Path("/nonexistent/index"),
                    "sha-current",
                    result_dir,
                    0.1,
                    0.1,
                    "direct",
                    result_dir,
                )


if __name__ == "__main__":
    unittest.main()
