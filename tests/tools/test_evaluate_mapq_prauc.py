#!/usr/bin/env python3
"""Tests for pair-level MAPQ PR-AUC evaluation."""

from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "evaluate-mapq-prauc.py"
SPEC = importlib.util.spec_from_file_location("bsbit_evaluate_mapq_prauc", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
PRAUC = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = PRAUC
SPEC.loader.exec_module(PRAUC)


class MapqPraucTests(unittest.TestCase):
    def test_ties_enter_together_and_missing_pairs_reduce_recall(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ledger = root / "ledger.tsv"
            with ledger.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("id", "state", "minimum_mapq", "proper"),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(
                    (
                        {"id": 0, "state": "exact", "minimum_mapq": 40, "proper": 1},
                        {
                            "id": 1,
                            "state": "wrong_axis",
                            "minimum_mapq": 40,
                            "proper": 1,
                        },
                        {"id": 2, "state": "exact", "minimum_mapq": 20, "proper": 1},
                        {"id": 3, "state": "absent", "minimum_mapq": -1, "proper": 0},
                        {"id": 4, "state": "near", "minimum_mapq": 0, "proper": 1},
                    )
                )
            args = PRAUC.parse_args(
                [
                    "--ledger",
                    str(ledger),
                    "--truth-pairs",
                    "5",
                    "--output-dir",
                    str(root / "result"),
                ]
            )
            summary = PRAUC.run(args)

            exact = summary["contracts"]["exact_origin"]
            locus = summary["contracts"]["locus_within_5bp"]
            self.assertAlmostEqual(exact["average_precision_step_prauc"], 7 / 30)
            self.assertAlmostEqual(exact["maximum_recall"], 2 / 5)
            self.assertAlmostEqual(locus["average_precision_step_prauc"], 23 / 60)
            self.assertAlmostEqual(locus["maximum_recall"], 3 / 5)
            self.assertEqual(
                exact["operating_points"]["q40"]["selected_pairs"], 2
            )
            self.assertAlmostEqual(
                exact["operating_points"]["q40"]["precision"], 1 / 2
            )
            self.assertAlmostEqual(
                exact["operating_points"]["q40"]["f1"], 2 / 7
            )
            self.assertAlmostEqual(
                locus["operating_points"]["q0"]["f1"], 2 / 3
            )

            with (root / "result" / "pr-curve.tsv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                rows = csv.DictReader(handle, delimiter="\t")
                self.assertIn("f1", rows.fieldnames or ())

    def test_truth_pair_count_must_match_ledger(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ledger = root / "ledger.tsv"
            ledger.write_text(
                "id\tstate\tminimum_mapq\tproper\n0\texact\t40\t1\n",
                encoding="utf-8",
            )
            with self.assertRaises(SystemExit):
                PRAUC.load_buckets(ledger, truth_pairs=2)


if __name__ == "__main__":
    unittest.main()
