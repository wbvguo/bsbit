#!/usr/bin/env python3
"""Generate deterministic directional WGBS pairs with known terminal contamination."""

from __future__ import annotations

import argparse
import random
from pathlib import Path


READ_BASES = 151
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def convert_top(sequence: str) -> str:
    return sequence.replace("C", "T")


def contaminated(sequence: str, clip: int, adapter: bool, rng: random.Random) -> str:
    if clip == 0:
        return sequence
    tail = (ADAPTER * ((clip + len(ADAPTER) - 1) // len(ADAPTER)))[:clip]
    if not adapter:
        tail = "".join(rng.choice("ACGT") for _ in range(clip))
    return sequence[:-clip] + tail


def top_strand_mismatch(*reference_bases: str) -> str:
    """Return a query base that is zero-cost against none of the references."""
    allowed = set(reference_bases)
    if "C" in reference_bases:
        allowed.add("T")
    return next(base for base in "ACGT" if base not in allowed)


def guaranteed_contaminant(
    reference_sequence: str,
    reverse: bool,
    alternate_reference_sequence: str | None = None,
) -> str:
    alternates = alternate_reference_sequence or reference_sequence
    if len(reference_sequence) != len(alternates):
        raise ValueError("reference alternatives must have equal lengths")
    oriented = "".join(
        top_strand_mismatch(primary, alternate)
        for primary, alternate in zip(reference_sequence, alternates, strict=True)
    )
    return reverse_complement(oriented) if reverse else oriented


def replace_terminals(
    sequence: str, contaminant: str, five_prime: int, three_prime: int
) -> str:
    if five_prime + three_prime >= len(sequence):
        raise ValueError("terminal contamination must retain at least one base")
    end = len(sequence) - three_prime if three_prime else len(sequence)
    return contaminant[:five_prime] + sequence[five_prime:end] + contaminant[end:]


def write_record(handle, name: str, sequence: str, quality: str) -> None:
    handle.write(f"@{name}\n{sequence}\n+\n{quality}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--pairs", type=int, default=20_000)
    parser.add_argument("--ambiguous-pairs", type=int, default=500)
    parser.add_argument("--ambiguous-clip", type=int, default=6)
    parser.add_argument("--reference-bases", type=int, default=2_000_000)
    parser.add_argument("--seed", type=int, default=0xB5B17)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.pairs <= 0:
        raise SystemExit("--pairs must be positive")
    if args.ambiguous_pairs < 0:
        raise SystemExit("--ambiguous-pairs must be nonnegative")
    if not 1 <= args.ambiguous_clip <= 30:
        raise SystemExit("--ambiguous-clip must be in 1..=30")
    if args.reference_bases < READ_BASES * 3:
        raise SystemExit("--reference-bases is too small")
    args.output.mkdir(parents=True, exist_ok=False)

    rng = random.Random(args.seed)
    reference_bases = [rng.choice("ACGT") for _ in range(args.reference_bases)]
    ambiguous_specs: list[tuple[int, int, int]] = []
    ambiguous_stride = 500
    ambiguous_fragment = 300
    ambiguous_clip = args.ambiguous_clip
    decoy_base = args.reference_bases // 2
    if decoy_base + args.ambiguous_pairs * ambiguous_stride + ambiguous_fragment >= len(
        reference_bases
    ):
        raise SystemExit("reference is too small for --ambiguous-pairs")
    for ordinal in range(args.ambiguous_pairs):
        truth_start = ordinal * ambiguous_stride
        decoy_start = decoy_base + ordinal * ambiguous_stride
        retained = READ_BASES - ambiguous_clip
        reference_bases[decoy_start : decoy_start + retained] = reference_bases[
            truth_start : truth_start + retained
        ]
        right_end = truth_start + ambiguous_fragment
        decoy_right_end = decoy_start + ambiguous_fragment
        reference_bases[decoy_right_end - retained : decoy_right_end] = reference_bases[
            right_end - retained : right_end
        ]
        ambiguous_specs.append((truth_start, decoy_start, ambiguous_fragment))
    reference = "".join(reference_bases)
    with (args.output / "reference.fa").open("w", encoding="ascii", newline="\n") as fasta:
        fasta.write(">truth\n")
        for start in range(0, len(reference), 80):
            fasta.write(reference[start : start + 80] + "\n")

    paths = {
        key: (args.output / name).open("w", encoding="ascii", newline="\n")
        for key, name in {
            "clean1": "clean_R1.fastq",
            "clean2": "clean_R2.fastq",
            "dirty1": "contaminated_R1.fastq",
            "dirty2": "contaminated_R2.fastq",
            "five1": "five_prime_R1.fastq",
            "five2": "five_prime_R2.fastq",
            "both1": "two_ended_R1.fastq",
            "both2": "two_ended_R2.fastq",
            "truth": "truth.tsv",
            "terminal_truth": "terminal_truth.tsv",
            "ambiguous_clean1": "ambiguous_clean_R1.fastq",
            "ambiguous_clean2": "ambiguous_clean_R2.fastq",
            "ambiguous_dirty1": "ambiguous_contaminated_R1.fastq",
            "ambiguous_dirty2": "ambiguous_contaminated_R2.fastq",
            "ambiguous_truth": "ambiguous_truth.tsv",
        }.items()
    }
    try:
        paths["truth"].write(
            "ordinal\tqname\tfragment_start\tfragment_length\t"
            "r1_start\tr2_start\tr1_clip\tr2_clip\n"
        )
        paths["ambiguous_truth"].write(
            "ordinal\tqname\ttruth_start\tdecoy_start\tfragment_length\tclip\n"
        )
        paths["terminal_truth"].write(
            "ordinal\tqname\tfive_r1\tthree_r1\tfive_r2\tthree_r2\n"
        )
        maximum_start = len(reference) - 500
        for ordinal in range(args.pairs):
            fragment_length = 260 + rng.randrange(161)
            fragment_start = rng.randrange(maximum_start - fragment_length)
            r1_reference = reference[fragment_start : fragment_start + READ_BASES]
            r2_start = fragment_start + fragment_length - READ_BASES
            r2_reference = reference[r2_start : r2_start + READ_BASES]
            read1 = convert_top(r1_reference)
            read2 = reverse_complement(convert_top(r2_reference))
            category = ordinal % 4
            r1_clip = 0 if category in (0, 2) else (30 if category == 1 else 25)
            r2_clip = 0 if category in (0, 1) else (30 if category == 2 else 25)
            dirty1 = contaminated(read1, r1_clip, True, rng)
            dirty2 = contaminated(read2, r2_clip, False, rng)
            contaminant1 = guaranteed_contaminant(r1_reference, reverse=False)
            contaminant2 = guaranteed_contaminant(r2_reference, reverse=True)
            five_r1 = 1 + ordinal % 6
            five_r2 = 1 + (ordinal // 6) % 6
            three_r1 = 1 + (ordinal // 36) % 6
            three_r2 = 1 + (ordinal // 216) % 6
            five_only1 = replace_terminals(read1, contaminant1, five_r1, 0)
            five_only2 = replace_terminals(read2, contaminant2, five_r2, 0)
            both1 = replace_terminals(read1, contaminant1, five_r1, three_r1)
            both2 = replace_terminals(read2, contaminant2, five_r2, three_r2)
            quality1 = "I" * READ_BASES
            quality2 = "I" * (READ_BASES - r2_clip) + "!" * r2_clip
            name = f"truth{ordinal}"
            write_record(paths["clean1"], name, read1, "I" * READ_BASES)
            write_record(paths["clean2"], name, read2, "I" * READ_BASES)
            write_record(paths["dirty1"], name, dirty1, quality1)
            write_record(paths["dirty2"], name, dirty2, quality2)
            write_record(paths["five1"], name, five_only1, "I" * READ_BASES)
            write_record(paths["five2"], name, five_only2, "I" * READ_BASES)
            write_record(paths["both1"], name, both1, "I" * READ_BASES)
            write_record(paths["both2"], name, both2, "I" * READ_BASES)
            paths["truth"].write(
                f"{ordinal}\t{name}\t{fragment_start}\t{fragment_length}\t"
                f"{fragment_start}\t{r2_start}\t{r1_clip}\t{r2_clip}\n"
            )
            paths["terminal_truth"].write(
                f"{ordinal}\t{name}\t{five_r1}\t{three_r1}\t"
                f"{five_r2}\t{three_r2}\n"
            )

        for ordinal, (truth_start, decoy_start, fragment_length) in enumerate(
            ambiguous_specs
        ):
            r2_start = truth_start + fragment_length - READ_BASES
            decoy_r2_start = decoy_start + fragment_length - READ_BASES
            r1_reference = reference[truth_start : truth_start + READ_BASES]
            decoy_r1_reference = reference[decoy_start : decoy_start + READ_BASES]
            r2_reference = reference[r2_start : r2_start + READ_BASES]
            decoy_r2_reference = reference[
                decoy_r2_start : decoy_r2_start + READ_BASES
            ]
            read1 = convert_top(r1_reference)
            read2 = reverse_complement(convert_top(r2_reference))
            dirty1 = replace_terminals(
                read1,
                guaranteed_contaminant(
                    r1_reference,
                    reverse=False,
                    alternate_reference_sequence=decoy_r1_reference,
                ),
                0,
                ambiguous_clip,
            )
            dirty2 = replace_terminals(
                read2,
                guaranteed_contaminant(
                    r2_reference,
                    reverse=True,
                    alternate_reference_sequence=decoy_r2_reference,
                ),
                0,
                ambiguous_clip,
            )
            quality = "I" * (READ_BASES - ambiguous_clip) + "!" * ambiguous_clip
            name = f"ambiguous{ordinal}"
            write_record(paths["ambiguous_clean1"], name, read1, "I" * READ_BASES)
            write_record(paths["ambiguous_clean2"], name, read2, "I" * READ_BASES)
            write_record(paths["ambiguous_dirty1"], name, dirty1, quality)
            write_record(paths["ambiguous_dirty2"], name, dirty2, quality)
            paths["ambiguous_truth"].write(
                f"{ordinal}\t{name}\t{truth_start}\t{decoy_start}\t"
                f"{fragment_length}\t{ambiguous_clip}\n"
            )
    finally:
        for handle in paths.values():
            handle.close()


if __name__ == "__main__":
    main()
