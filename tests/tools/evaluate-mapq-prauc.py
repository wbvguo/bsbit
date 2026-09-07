#!/usr/bin/env python3
"""Compute pair-level mapping PR-AUC from a simulator-truth audit ledger.

The score is the minimum MAPQ of the two primary mates.  Equal MAPQ values
enter the curve as one group, so the result cannot depend on an arbitrary
ordering within a tie.  Recall is divided by every input truth pair rather
than only by reported mappings; missing and non-proper pairs therefore earn
no recall and no area.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


SCHEMA = "bsbit-pair-mapq-prauc"
MAX_MAPQ = 60
CONTRACTS = {
    "exact_origin": frozenset({"exact"}),
    "locus_within_5bp": frozenset({"exact", "near"}),
}
PROPER_STATES = frozenset({"exact", "near", "same_axis", "wrong_axis"})
REQUIRED_COLUMNS = frozenset({"id", "state", "minimum_mapq", "proper"})


def fail(message: str) -> None:
    raise SystemExit(f"error: {message}")


def parse_positive(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected an integer, found {value!r}") from None
    if parsed <= 0:
        raise argparse.ArgumentTypeError(f"expected a positive integer, found {parsed}")
    return parsed


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--truth-pairs", required=True, type=parse_positive)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args(argv)


@dataclass(frozen=True)
class ScoreBucket:
    reported: int
    states: dict[str, int]


@dataclass(frozen=True)
class CurvePoint:
    threshold: int
    selected: int
    correct: int
    incorrect: int
    precision: float
    recall: float
    f1: float


def load_buckets(path: Path, truth_pairs: int) -> tuple[dict[int, ScoreBucket], int]:
    counts: dict[int, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    reported_by_score: dict[int, int] = defaultdict(int)
    row_count = 0
    proper_count = 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if rows.fieldnames is None:
            fail(f"ledger {path} has no header")
        missing = REQUIRED_COLUMNS.difference(rows.fieldnames)
        if missing:
            fail(f"ledger {path} is missing columns: {', '.join(sorted(missing))}")
        for row_count, row in enumerate(rows, start=1):
            context = f"{path}:{row_count + 1}"
            try:
                pair_id = int(row["id"])
                proper = int(row["proper"])
                mapq = int(row["minimum_mapq"])
            except ValueError:
                fail(f"{context}: id, proper, and minimum_mapq must be integers")
            if pair_id != row_count - 1:
                fail(
                    f"{context}: expected ordered pair id {row_count - 1}, found {pair_id}"
                )
            if proper not in (0, 1):
                fail(f"{context}: proper must be zero or one, found {proper}")
            if not proper:
                if mapq != -1:
                    fail(f"{context}: non-proper pair must have minimum_mapq -1")
                continue
            state = row["state"]
            if state not in PROPER_STATES:
                fail(f"{context}: proper pair has unsupported state {state!r}")
            if not 0 <= mapq <= MAX_MAPQ:
                fail(f"{context}: proper pair MAPQ must be in 0..={MAX_MAPQ}, found {mapq}")
            proper_count += 1
            reported_by_score[mapq] += 1
            counts[mapq][state] += 1
    if row_count != truth_pairs:
        fail(f"ledger contains {row_count} pairs but --truth-pairs is {truth_pairs}")
    buckets = {
        score: ScoreBucket(reported=reported_by_score[score], states=dict(states))
        for score, states in counts.items()
    }
    return buckets, proper_count


def build_curve(
    buckets: dict[int, ScoreBucket], correct_states: frozenset[str], truth_pairs: int
) -> tuple[list[CurvePoint], float, float]:
    selected = 0
    correct = 0
    previous_precision = 1.0
    previous_recall = 0.0
    average_precision = 0.0
    trapezoidal_auc = 0.0
    points = [CurvePoint(MAX_MAPQ + 1, 0, 0, 0, 1.0, 0.0, 0.0)]
    for threshold in range(MAX_MAPQ, -1, -1):
        bucket = buckets.get(threshold)
        if bucket is not None:
            selected += bucket.reported
            correct += sum(bucket.states.get(state, 0) for state in correct_states)
        precision = correct / selected if selected else 1.0
        recall = correct / truth_pairs
        f1 = (
            2.0 * correct / (truth_pairs + selected)
            if truth_pairs + selected
            else 0.0
        )
        recall_increment = recall - previous_recall
        average_precision += recall_increment * precision
        trapezoidal_auc += recall_increment * (previous_precision + precision) / 2.0
        points.append(
            CurvePoint(
                threshold=threshold,
                selected=selected,
                correct=correct,
                incorrect=selected - correct,
                precision=precision,
                recall=recall,
                f1=f1,
            )
        )
        previous_precision = precision
        previous_recall = recall
    return points, average_precision, trapezoidal_auc


def finite(value: float, context: str) -> float:
    if not math.isfinite(value):
        fail(f"non-finite result for {context}")
    return value


def write_curve(path: Path, named_points: Iterable[tuple[str, CurvePoint]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "contract",
                "mapq_threshold",
                "selected_pairs",
                "correct_pairs",
                "incorrect_pairs",
                "precision",
                "recall",
                "f1",
            )
        )
        for contract, point in named_points:
            writer.writerow(
                (
                    contract,
                    point.threshold,
                    point.selected,
                    point.correct,
                    point.incorrect,
                    f"{point.precision:.12g}",
                    f"{point.recall:.12g}",
                    f"{point.f1:.12g}",
                )
            )


def run(args: argparse.Namespace) -> dict[str, object]:
    ledger = args.ledger.resolve(strict=True)
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        fail(f"refusing to overwrite output directory {output_dir}")
    output_dir.mkdir(parents=True)
    buckets, proper_count = load_buckets(ledger, args.truth_pairs)
    summaries: dict[str, object] = {}
    named_points: list[tuple[str, CurvePoint]] = []
    for contract, correct_states in CONTRACTS.items():
        points, average_precision, trapezoidal_auc = build_curve(
            buckets, correct_states, args.truth_pairs
        )
        named_points.extend((contract, point) for point in points)
        final = points[-1]
        by_threshold = {point.threshold: point for point in points}
        summaries[contract] = {
            "correct_reported_pairs": final.correct,
            "incorrect_reported_pairs": final.incorrect,
            "maximum_recall": finite(final.recall, f"{contract} maximum recall"),
            "all_reported_precision": finite(
                final.precision, f"{contract} all-reported precision"
            ),
            "average_precision_step_prauc": finite(
                average_precision, f"{contract} average precision"
            ),
            "trapezoidal_prauc": finite(
                trapezoidal_auc, f"{contract} trapezoidal PR-AUC"
            ),
            "operating_points": {
                f"q{threshold}": {
                    "selected_pairs": by_threshold[threshold].selected,
                    "correct_pairs": by_threshold[threshold].correct,
                    "incorrect_pairs": by_threshold[threshold].incorrect,
                    "precision": by_threshold[threshold].precision,
                    "recall": by_threshold[threshold].recall,
                    "f1": by_threshold[threshold].f1,
                }
                for threshold in (0, 10, 20, 30, 40)
            },
        }
    summary: dict[str, object] = {
        "schema": SCHEMA,
        "ledger": str(ledger),
        "truth_pairs": args.truth_pairs,
        "reported_proper_pairs": proper_count,
        "score": "minimum MAPQ of the two primary mates",
        "tie_policy": "all pairs with equal integer MAPQ enter simultaneously",
        "recall_denominator": "all input truth pairs",
        "primary_metric": "average_precision_step_prauc",
        "contracts": summaries,
    }
    write_curve(output_dir / "pr-curve.tsv", named_points)
    with (output_dir / "summary.json").open("x", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def main(argv: Sequence[str] | None = None) -> None:
    summary = run(parse_args(argv))
    print("contract\tstep_prauc\ttrapezoidal_prauc\tmax_recall\tall_precision")
    contracts = summary["contracts"]
    assert isinstance(contracts, dict)
    for contract in CONTRACTS:
        values = contracts[contract]
        assert isinstance(values, dict)
        print(
            contract,
            f'{values["average_precision_step_prauc"]:.12g}',
            f'{values["trapezoidal_prauc"]:.12g}',
            f'{values["maximum_recall"]:.12g}',
            f'{values["all_reported_precision"]:.12g}',
            sep="\t",
        )


if __name__ == "__main__":
    main()
