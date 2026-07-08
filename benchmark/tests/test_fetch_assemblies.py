from __future__ import annotations

import io
import unittest
from contextlib import redirect_stderr

from sparrowhawk_amr_benchmark.fetch_assemblies import print_final_missing_ids


class PrintFinalMissingIdsTests(unittest.TestCase):
    def test_prints_missing_assembly_ids_in_ordered_row_order(self) -> None:
        rows = [
            {"species": "Zeta", "assembly_id": "GCF_000003"},
            {"species": "Alpha", "assembly_id": "GCF_000002"},
            {"species": "Alpha", "assembly_id": "GCA_000001"},
        ]

        stderr = io.StringIO()
        with redirect_stderr(stderr):
            print_final_missing_ids(rows)

        self.assertEqual(
            stderr.getvalue(),
            "\n".join(
                [
                    "Assemblies still missing after all retries:",
                    "  GCA_000001",
                    "  GCF_000002",
                    "  GCF_000003",
                    "",
                ]
            ),
        )

    def test_prints_nothing_when_no_rows_are_missing(self) -> None:
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            print_final_missing_ids([])

        self.assertEqual(stderr.getvalue(), "")


if __name__ == "__main__":
    unittest.main()
