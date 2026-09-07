#!/usr/bin/env python3
"""Evaluate ordered single-end primary alignments against paired-read truth.

The truth BAM contributes primary read 1 only.  Candidate output must preserve
that input order, but may omit records.  Coordinates use the unclipped SAM
five-prime origin: the left unclipped edge on a forward record and the right
unclipped edge on a reverse record.  Equal integer MAPQ scores enter the
precision/recall curve together, and recall always uses every truth read.

Only the Python standard library and ``samtools view`` are required.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence


TOLERANCES = (0, 1, 5, 10)
PRIMARY_EXCLUDED_FLAGS = 0x100 | 0x800
READ1_FLAG = 0x40
UNMAPPED_FLAG = 0x4
REVERSE_FLAG = 0x10
CIGAR_TOKEN = re.compile(r"([0-9]+)([MIDNSHP=X])")
REFERENCE_CONSUMING = frozenset("MDN=X")
CLIPPING = frozenset("SH")


def fail(message: str) -> None:
    raise SystemExit(f"error: {message}")


def parse_nonnegative(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected an integer, found {value!r}") from None
    if parsed < 0:
        raise argparse.ArgumentTypeError(f"expected a nonnegative integer, found {parsed}")
    return parsed


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth-bam", required=True, type=Path)
    parser.add_argument("--candidate-bam", required=True, type=Path)
    parser.add_argument("--aligner", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--maximum-mapq", type=parse_nonnegative, default=60)
    return parser.parse_args(argv)


@dataclass(frozen=True)
class SamRecord:
    qname: str
    flag: int
    rname: str
    position: int
    mapq: int
    cigar: str

    @classmethod
    def parse(cls, line: str) -> "SamRecord":
        fields = line.rstrip("\n").split("\t", 11)
        if len(fields) < 6:
            raise ValueError("SAM record has fewer than six fields")
        try:
            return cls(
                qname=fields[0],
                flag=int(fields[1]),
                rname=fields[2],
                position=int(fields[3]),
                mapq=int(fields[4]),
                cigar=fields[5],
            )
        except ValueError as error:
            raise ValueError("SAM FLAG, POS, and MAPQ must be integers") from error

    @property
    def primary(self) -> bool:
        return self.flag & PRIMARY_EXCLUDED_FLAGS == 0

    @property
    def read1(self) -> bool:
        return self.flag & READ1_FLAG != 0

    @property
    def mapped(self) -> bool:
        return self.flag & UNMAPPED_FLAG == 0

    @property
    def reverse(self) -> bool:
        return self.flag & REVERSE_FLAG != 0

    def five_prime_origin(self) -> int:
        if not self.mapped or self.position <= 0 or self.cigar == "*":
            raise ValueError(f"mapped record {self.qname!r} has no coordinate-bearing CIGAR")
        operations = parse_cigar(self.cigar)
        reference_length = sum(
            length for length, operation in operations if operation in REFERENCE_CONSUMING
        )
        if reference_length == 0:
            raise ValueError(f"mapped record {self.qname!r} consumes no reference bases")
        if self.reverse:
            trailing_clip = sum_trailing_clip(operations)
            return self.position + reference_length - 1 + trailing_clip
        leading_clip = sum_leading_clip(operations)
        return self.position - leading_clip


def parse_cigar(cigar: str) -> tuple[tuple[int, str], ...]:
    operations = tuple(
        (int(length), operation) for length, operation in CIGAR_TOKEN.findall(cigar)
    )
    if not operations or "".join(f"{length}{operation}" for length, operation in operations) != cigar:
        raise ValueError(f"invalid CIGAR {cigar!r}")
    if any(length <= 0 for length, _ in operations):
        raise ValueError(f"CIGAR contains a zero-length operation: {cigar!r}")
    return operations


def sum_leading_clip(operations: Sequence[tuple[int, str]]) -> int:
    total = 0
    for length, operation in operations:
        if operation not in CLIPPING:
            break
        total += length
    return total


def sum_trailing_clip(operations: Sequence[tuple[int, str]]) -> int:
    total = 0
    for length, operation in reversed(operations):
        if operation not in CLIPPING:
            break
        total += length
    return total


def normalized_qname(qname: str) -> str:
    return qname[:-2] if qname.endswith(("/1", "/2")) else qname


def samtools_records(command: str, path: Path) -> Iterator[SamRecord]:
    resolved = path.resolve(strict=True)
    process = subprocess.Popen(
        [command, "view", str(resolved)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )
    assert process.stdout is not None and process.stderr is not None
    try:
        for line_number, line in enumerate(process.stdout, start=1):
            try:
                yield SamRecord.parse(line)
            except ValueError as error:
                fail(f"{resolved}: SAM record {line_number}: {error}")
    finally:
        process.stdout.close()
    stderr = process.stderr.read()
    status = process.wait()
    if status:
        fail(f"{command} view {resolved} failed with status {status}: {stderr.strip()}")


def truth_read1(records: Iterable[SamRecord]) -> Iterator[SamRecord]:
    for record in records:
        if record.primary and record.read1:
            if not record.mapped:
                fail(f"truth read {record.qname!r} is unmapped")
            yield record


def candidate_primary(records: Iterable[SamRecord]) -> Iterator[SamRecord]:
    yield from (record for record in records if record.primary)


def average_precision(
    reported_by_score: Counter[int],
    correct_by_score: Counter[int],
    input_reads: int,
    maximum_mapq: int,
) -> float:
    selected = 0
    correct = 0
    previous_recall = 0.0
    area = 0.0
    for threshold in range(maximum_mapq, -1, -1):
        selected += reported_by_score[threshold]
        correct += correct_by_score[threshold]
        precision = correct / selected if selected else 1.0
        recall = correct / input_reads if input_reads else 0.0
        area += (recall - previous_recall) * precision
        previous_recall = recall
    return area


def evaluate(
    truth: Iterable[SamRecord],
    candidate: Iterable[SamRecord],
    aligner: str,
    maximum_mapq: int,
    progress: bool = False,
) -> dict[str, object]:
    truth_iterator = iter(truth)
    candidate_iterator = iter(candidate)
    pending = next(candidate_iterator, None)
    input_reads = 0
    primary_records = 0
    mapped = 0
    omitted = 0
    same_axis = 0
    same_axis_over_10bp = 0
    wrong_contig = 0
    wrong_orientation = 0
    unknown_mapq = 0
    above_maximum_mapq = 0
    mapped_mapq_ge_1 = 0
    raw_mapq_counts: Counter[int] = Counter()
    reported_by_score: Counter[int] = Counter()
    correct_by_score = {tolerance: Counter() for tolerance in TOLERANCES}
    correct_counts = Counter()

    for expected in truth_iterator:
        input_reads += 1
        if progress and input_reads % 1_000_000 == 0:
            print(f"truth R1 records: {input_reads:,}", file=sys.stderr)
        expected_name = normalized_qname(expected.qname)
        if pending is None or normalized_qname(pending.qname) != expected_name:
            omitted += 1
            continue
        observed = pending
        pending = next(candidate_iterator, None)
        primary_records += 1
        if progress and primary_records % 1_000_000 == 0:
            print(f"candidate primary records: {primary_records:,}", file=sys.stderr)
        raw_mapq_counts[observed.mapq] += 1
        if not observed.mapped:
            continue
        mapped += 1
        if observed.mapq == 255:
            unknown_mapq += 1
            effective_mapq = 0
        else:
            if observed.mapq > maximum_mapq:
                above_maximum_mapq += 1
            effective_mapq = min(maximum_mapq, max(0, observed.mapq))
        mapped_mapq_ge_1 += int(effective_mapq >= 1)
        reported_by_score[effective_mapq] += 1
        if observed.rname != expected.rname:
            wrong_contig += 1
            continue
        if observed.reverse != expected.reverse:
            wrong_orientation += 1
            continue
        same_axis += 1
        displacement = abs(observed.five_prime_origin() - expected.five_prime_origin())
        same_axis_over_10bp += int(displacement > 10)
        for tolerance in TOLERANCES:
            if displacement <= tolerance:
                correct_counts[tolerance] += 1
                correct_by_score[tolerance][effective_mapq] += 1

    if pending is not None:
        fail(
            "candidate primary order diverges from truth at "
            f"{pending.qname!r}; candidate output must preserve input order"
        )
    trailing = next(candidate_iterator, None)
    if trailing is not None:
        fail(f"candidate contains an unexpected trailing primary record {trailing.qname!r}")
    if input_reads == 0:
        fail("truth BAM contains no primary read-1 records")

    summaries: dict[str, object] = {}
    for tolerance in TOLERANCES:
        correct = correct_counts[tolerance]
        errors = mapped - correct
        precision = correct / mapped if mapped else 1.0
        recall = correct / input_reads
        f1 = (
            2.0 * precision * recall / (precision + recall)
            if precision + recall
            else 0.0
        )
        summaries[str(tolerance)] = {
            "average_precision_step_prauc": average_precision(
                reported_by_score,
                correct_by_score[tolerance],
                input_reads,
                maximum_mapq,
            ),
            "correct": correct,
            "errors": errors,
            "f1": f1,
            "precision": precision,
            "recall": recall,
            "reported": mapped,
        }
    return {
        "accuracy_by_tolerance_bp": summaries,
        "aligner": aligner,
        "axis_policy_requested": "sam",
        "axis_policy_resolved": "sam",
        "counts": {
            "mapped_effective_mapq_ge_1": mapped_mapq_ge_1,
            "mapped_primary_records": mapped,
            "mapq_above_max_mapped_records": above_maximum_mapq,
            "omitted_input_qnames": omitted,
            "primary_records": primary_records,
            "records": primary_records,
            "same_axis": same_axis,
            "same_axis_over_10bp": same_axis_over_10bp,
            "unknown_mapq_mapped_records": unknown_mapq,
            **{f"within_{tolerance}bp_mapped": correct_counts[tolerance] for tolerance in TOLERANCES},
            "wrong_contig": wrong_contig,
            "wrong_orientation": wrong_orientation,
        },
        "input_reads": input_reads,
        "maximum_mapq": maximum_mapq,
        "raw_mapq_counts": {
            str(score): count for score, count in sorted(raw_mapq_counts.items())
        },
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        fail(f"refusing to overwrite output directory {output_dir}")
    summary = evaluate(
        truth_read1(samtools_records(args.samtools, args.truth_bam)),
        candidate_primary(samtools_records(args.samtools, args.candidate_bam)),
        args.aligner,
        args.maximum_mapq,
        progress=True,
    )
    output_dir.mkdir(parents=True)
    with (output_dir / "summary.json").open("x", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def main(argv: Sequence[str] | None = None) -> None:
    summary = run(parse_args(argv))
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
