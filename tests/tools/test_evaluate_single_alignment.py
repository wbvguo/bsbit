#!/usr/bin/env python3
"""Tests for ordered single-end truth evaluation."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "tests" / "tools" / "evaluate-single-alignment.py"
SPEC = importlib.util.spec_from_file_location("bsbit_evaluate_single_alignment", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
EVALUATOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = EVALUATOR
SPEC.loader.exec_module(EVALUATOR)


def record(
    qname: str,
    flag: int,
    rname: str,
    position: int,
    mapq: int,
    cigar: str = "100M",
):
    return EVALUATOR.SamRecord(qname, flag, rname, position, mapq, cigar)


class SingleAlignmentEvaluationTests(unittest.TestCase):
    def test_unclipped_five_prime_origin_uses_strand_specific_edge(self) -> None:
        forward = record("forward", 0, "chr1", 101, 20, "3H5S95M")
        reverse = record("reverse", 0x10, "chr1", 101, 20, "95M5S3H")
        self.assertEqual(forward.five_prime_origin(), 93)
        self.assertEqual(reverse.five_prime_origin(), 203)

    def test_tied_mapq_curve_and_axis_categories_use_all_truth_reads(self) -> None:
        truth = [
            record("r0", 0x40, "chr1", 100, 60),
            record("r1", 0x40 | 0x10, "chr1", 200, 60),
            record("r2", 0x40, "chr1", 300, 60),
            record("r3", 0x40, "chr1", 400, 60),
            record("r4", 0x40, "chr2", 500, 60),
        ]
        candidate = [
            record("r0/1", 0, "chr1", 100, 40),
            record("r1/1", 0x10, "chr1", 196, 20),
            record("r2/1", 0x10, "chr1", 300, 30),
            record("r3/1", 0x4, "*", 0, 0, "*"),
            record("r4/1", 0, "chr1", 500, 10),
        ]

        summary = EVALUATOR.evaluate(truth, candidate, "fixture", 60)
        counts = summary["counts"]
        self.assertEqual(counts["mapped_primary_records"], 4)
        self.assertEqual(counts["same_axis"], 2)
        self.assertEqual(counts["wrong_orientation"], 1)
        self.assertEqual(counts["wrong_contig"], 1)
        self.assertEqual(counts["within_0bp_mapped"], 1)
        self.assertEqual(counts["within_5bp_mapped"], 2)
        within_five = summary["accuracy_by_tolerance_bp"]["5"]
        self.assertAlmostEqual(within_five["average_precision_step_prauc"], 1 / 3)
        self.assertAlmostEqual(within_five["precision"], 1 / 2)
        self.assertAlmostEqual(within_five["recall"], 2 / 5)
        self.assertAlmostEqual(within_five["f1"], 4 / 9)

    def test_ordered_candidate_omission_remains_in_recall_denominator(self) -> None:
        truth = [
            record("missing", 0x40, "chr1", 10, 60),
            record("reported", 0x40, "chr1", 20, 60),
        ]
        candidate = [record("reported/1", 0, "chr1", 20, 40)]
        summary = EVALUATOR.evaluate(truth, candidate, "fixture", 60)
        self.assertEqual(summary["input_reads"], 2)
        self.assertEqual(summary["counts"]["omitted_input_qnames"], 1)
        self.assertEqual(summary["counts"]["primary_records"], 1)
        self.assertAlmostEqual(
            summary["accuracy_by_tolerance_bp"]["0"]["recall"], 1 / 2
        )


if __name__ == "__main__":
    unittest.main()
