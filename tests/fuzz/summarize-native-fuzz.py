#!/usr/bin/env python3
"""Validate and summarize bounded native libFuzzer logs."""

from __future__ import annotations

import argparse
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


SCHEMA = "bsbit-native-coverage-fuzz"
TARGETS = ("reader", "record", "header")
DONE = re.compile(
    r"^#(?P<executions>\d+)\s+DONE\s+cov:\s+(?P<coverage>\d+)\s+"
    r"ft:\s+(?P<features>\d+)\s+corp:\s+(?P<files>\d+)/"
    r"(?P<bytes>\d+)(?P<suffix>[KMG]?)b\b.*\srss:\s+(?P<rss>\d+)Mb$"
)
ELAPSED = re.compile(
    r"^Done (?P<executions>\d+) runs in (?P<seconds>\d+) second\(s\)$"
)
EXECUTED = re.compile(r"^stat::number_of_executed_units:\s+(?P<value>\d+)$")
PEAK_RSS = re.compile(r"^stat::peak_rss_mb:\s+(?P<value>\d+)$")
FAILURE_MARKERS = (
    "ERROR: AddressSanitizer",
    "ERROR: LeakSanitizer",
    "ERROR: libFuzzer",
    "SUMMARY: UndefinedBehaviorSanitizer",
    "LeakSanitizer has encountered a fatal error",
    "native fuzz requirement failed",
    "runtime error:",
)


class ValidationError(RuntimeError):
    """One native fuzz log is incomplete or internally inconsistent."""


@dataclass(frozen=True)
class Summary:
    """Validated final statistics for one target."""

    target: str
    executions: int
    coverage: int
    features: int
    corpus_files: int
    corpus_bytes: int
    elapsed_seconds: int
    peak_rss_mb: int


def one_match(
    pattern: re.Pattern[str], lines: list[str], label: str
) -> re.Match[str]:
    """Return one required regex match."""

    matches = [match for line in lines if (match := pattern.match(line))]
    if len(matches) != 1:
        raise ValidationError(f"expected exactly one {label}, found {len(matches)}")
    return matches[0]


def corpus_bytes(value: str, suffix: str) -> int:
    """Expand libFuzzer's binary corpus-size suffix."""

    multiplier = {"": 1, "K": 1_024, "M": 1_048_576, "G": 1_073_741_824}[
        suffix
    ]
    return int(value) * multiplier


def parse_log(path: Path, requested_seconds: int) -> Summary:
    """Validate one native libFuzzer log."""

    target = path.stem
    if target not in TARGETS:
        raise ValidationError(f"unexpected native fuzz target log: {target}")
    text = path.read_text(encoding="utf-8", errors="replace")
    for marker in FAILURE_MARKERS:
        if marker in text:
            raise ValidationError(f"{target} log contains failure marker: {marker}")
    lines = text.splitlines()
    done = one_match(DONE, lines, "DONE line")
    elapsed = one_match(ELAPSED, lines, "elapsed line")
    executed = one_match(EXECUTED, lines, "executed-units statistic")
    peak = one_match(PEAK_RSS, lines, "peak-RSS statistic")
    execution_values = {
        int(done.group("executions")),
        int(elapsed.group("executions")),
        int(executed.group("value")),
    }
    if len(execution_values) != 1:
        raise ValidationError(f"{target} execution counts disagree: {execution_values}")
    rss_values = {int(done.group("rss")), int(peak.group("value"))}
    if len(rss_values) != 1:
        raise ValidationError(f"{target} peak RSS values disagree: {rss_values}")
    elapsed_seconds = int(elapsed.group("seconds"))
    if not requested_seconds <= elapsed_seconds <= requested_seconds + 5:
        raise ValidationError(
            f"{target} elapsed seconds {elapsed_seconds} do not prove "
            f"requested budget {requested_seconds}"
        )
    values = {
        "executions": execution_values.pop(),
        "coverage": int(done.group("coverage")),
        "features": int(done.group("features")),
        "corpus_files": int(done.group("files")),
        "peak_rss_mb": rss_values.pop(),
    }
    for label, value in values.items():
        if value <= 0:
            raise ValidationError(f"{target} {label} must be positive, found {value}")
    return Summary(
        target=target,
        executions=values["executions"],
        coverage=values["coverage"],
        features=values["features"],
        corpus_files=values["corpus_files"],
        corpus_bytes=corpus_bytes(done.group("bytes"), done.group("suffix")),
        elapsed_seconds=elapsed_seconds,
        peak_rss_mb=values["peak_rss_mb"],
    )


def write_summary(
    logs: list[Path],
    output: Path,
    requested_seconds: int,
    seed: int,
    max_len: int,
) -> None:
    """Write one exact TSV row per declared target."""

    summaries = sorted(
        (parse_log(path, requested_seconds) for path in logs),
        key=lambda row: TARGETS.index(row.target),
    )
    if tuple(row.target for row in summaries) != TARGETS:
        raise ValidationError("summary requires each declared target exactly once")
    lines = [
        "\t".join(
            (
                "schema",
                "target",
                "requested_seconds",
                "executed_units",
                "coverage_edges",
                "features",
                "corpus_files",
                "corpus_bytes",
                "elapsed_seconds",
                "peak_rss_mb",
                "seed",
                "max_len",
            )
        )
    ]
    for row in summaries:
        lines.append(
            "\t".join(
                str(value)
                for value in (
                    SCHEMA,
                    row.target,
                    requested_seconds,
                    row.executions,
                    row.coverage,
                    row.features,
                    row.corpus_files,
                    row.corpus_bytes,
                    row.elapsed_seconds,
                    row.peak_rss_mb,
                    seed,
                    max_len,
                )
            )
        )
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def self_test() -> None:
    """Exercise positive parsing and required negative mutations."""

    template = """#123\tDONE   cov: 45 ft: 67 corp: 8/2Kb lim: 10 exec/s: 41 rss: 52Mb
Done 123 runs in 4 second(s)
stat::number_of_executed_units: 123
stat::peak_rss_mb:              52
"""
    with tempfile.TemporaryDirectory(prefix="bsbit-native-fuzz-summary-") as directory:
        root = Path(directory)
        logs = []
        for target in TARGETS:
            path = root / f"{target}.log"
            path.write_text(template, encoding="utf-8")
            logs.append(path)
        output = root / "summary.tsv"
        write_summary(logs, output, 3, 7, 4_096)
        rows = output.read_text(encoding="utf-8").splitlines()
        if len(rows) != 4 or "\treader\t" not in rows[1] or "\t2048\t" not in rows[1]:
            raise ValidationError("positive native summary self-test failed")

        logs[0].write_text(template.replace("Done 123", "Done 124"), encoding="utf-8")
        try:
            write_summary(logs, output, 3, 7, 4_096)
        except ValidationError:
            pass
        else:
            raise ValidationError("execution-count mutation was accepted")

        logs[0].write_text(template + "ERROR: AddressSanitizer\n", encoding="utf-8")
        try:
            write_summary(logs, output, 3, 7, 4_096)
        except ValidationError:
            pass
        else:
            raise ValidationError("sanitizer marker was accepted")

        logs[0].write_text(template.replace("in 4 second", "in 2 second"), encoding="utf-8")
        try:
            write_summary(logs, output, 3, 7, 4_096)
        except ValidationError:
            pass
        else:
            raise ValidationError("under-duration log was accepted")
    print("Native fuzz summary self-test passed")


def main() -> int:
    """Run the selected validator command."""

    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    summary = commands.add_parser("summary")
    summary.add_argument("--requested-seconds", type=int, required=True)
    summary.add_argument("--seed", type=int, required=True)
    summary.add_argument("--max-len", type=int, required=True)
    summary.add_argument("--output", type=Path, required=True)
    summary.add_argument("logs", nargs="+", type=Path)
    commands.add_parser("self-test")
    args = parser.parse_args()
    try:
        if args.command == "self-test":
            self_test()
        else:
            if args.requested_seconds < 1:
                raise ValidationError("requested seconds must be positive")
            write_summary(
                args.logs,
                args.output,
                args.requested_seconds,
                args.seed,
                args.max_len,
            )
    except (OSError, ValidationError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
