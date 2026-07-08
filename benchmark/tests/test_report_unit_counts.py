from __future__ import annotations

import unittest

from sparrowhawk_amr_benchmark.compare_amrfinder_batch import (
    report_map_context,
    report_unit_counts,
)


class ReportUnitCountsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.context = report_map_context(
            [
                {"type": "AMR", "hierarchy_node": "parent"},
                {"type": "AMR", "hierarchy_node": "child_a"},
                {"type": "AMR", "hierarchy_node": "child_b"},
                {"type": "AMR", "hierarchy_node": "child_c"},
                {"type": "AMR", "hierarchy_node": "unrelated_parent"},
            ],
            {
                "parent": {"parent_node_id": ""},
                "child_a": {"parent_node_id": "parent"},
                "child_b": {"parent_node_id": "parent"},
                "child_c": {"parent_node_id": "parent"},
                "unrelated_parent": {"parent_node_id": ""},
            },
            {"AMR"},
        )
        self.universe = {"parent", "child_a", "child_b", "child_c", "unrelated_parent"}

    def test_parent_detector_covers_multiple_truth_children_without_descendant_fp(self) -> None:
        counts = report_unit_counts({"parent"}, {"child_a", "child_b"}, self.universe, self.context)

        self.assertEqual(counts["tp"], 2)
        self.assertEqual(counts["fn"], 0)
        self.assertEqual(counts["fp"], 0)

    def test_parent_detector_covers_single_truth_child_without_unreported_descendant_fp(self) -> None:
        counts = report_unit_counts({"parent"}, {"child_a"}, self.universe, self.context)

        self.assertEqual(counts["tp"], 1)
        self.assertEqual(counts["fn"], 0)
        self.assertEqual(counts["fp"], 0)

    def test_unrelated_parent_is_detector_only_and_truth_is_missed(self) -> None:
        counts = report_unit_counts({"unrelated_parent"}, {"child_a"}, self.universe, self.context)

        self.assertEqual(counts["tp"], 0)
        self.assertEqual(counts["fn"], 1)
        self.assertEqual(counts["fp"], 1)


if __name__ == "__main__":
    unittest.main()
