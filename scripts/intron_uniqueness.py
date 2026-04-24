#!/usr/bin/env python3
"""Summarize how uniquely annotated intronic sequence maps to genes.

This script reads a GTF, derives gene-level introns from gene + exon features,
and reports how much intronic sequence is unique to one gene:

- on the same strand
- ignoring strand

The calculation mirrors the one used to interpret the nuclei RNA-seq dataset in
this repo. It uses gene spans, not transcript models, so an intronic base is
treated as ambiguous if it overlaps any other annotated gene body.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GTF = REPO_ROOT / "resources" / "reference" / "Mus_musculus.GRCm39.115.gtf.gz"
DEFAULT_OUTPUT = REPO_ROOT / "reports" / "intron_uniqueness_metrics.csv"

### COMMAND-LINE ARGUMENTS ####################################################


def parse_args() -> argparse.Namespace:
    """Parse the small set of inputs needed for this one-off annotation scan."""
    parser = argparse.ArgumentParser(
        description="Estimate how uniquely annotated intronic sequence maps to genes."
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        default=DEFAULT_GTF,
        help=f"GTF file to analyze. Default: {DEFAULT_GTF}",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"CSV file to write. Default: {DEFAULT_OUTPUT}",
    )
    return parser.parse_args()

### GTF PARSING ###############################################################


def parse_gtf_attributes(attribute_text: str) -> dict[str, str]:
    """Extract `key "value"` pairs from the ninth GTF column."""
    attributes: dict[str, str] = {}
    for part in attribute_text.strip().split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, value = part.split(" ", 1)
        attributes[key] = value.strip().strip('"')
    return attributes


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or abutting genomic intervals."""
    if not intervals:
        return []

    sorted_intervals = sorted(intervals)
    merged: list[list[int]] = [list(sorted_intervals[0])]
    for start, end in sorted_intervals[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def load_gene_and_exon_annotations(
    gtf_path: Path,
) -> tuple[dict[str, tuple[str, str, int, int]], dict[str, list[tuple[int, int]]]]:
    """Load gene spans and exon intervals keyed by gene ID."""
    genes: dict[str, tuple[str, str, int, int]] = {}
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)

    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            chrom, _source, feature, start, end, _score, strand, _frame, attrs = line.rstrip(
                "\n"
            ).split("\t")

            if feature not in {"gene", "exon"}:
                continue

            gene_id = parse_gtf_attributes(attrs).get("gene_id")
            if not gene_id:
                continue

            start_i = int(start)
            end_i = int(end)

            if feature == "gene":
                genes[gene_id] = (chrom, strand, start_i, end_i)
            else:
                exons[gene_id].append((start_i, end_i))

    return genes, exons

### INTRON DERIVATION #########################################################


def derive_introns(
    genes: dict[str, tuple[str, str, int, int]],
    exons: dict[str, list[tuple[int, int]]],
) -> tuple[
    dict[str, list[tuple[int, int, str, str]]],
    dict[str, list[tuple[int, int, str, str]]],
]:
    """Return introns and gene spans grouped by chromosome."""

    introns_by_chr: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)
    gene_spans_by_chr: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)

    for gene_id, (chrom, strand, gene_start, gene_end) in genes.items():
        gene_spans_by_chr[chrom].append((gene_start, gene_end, strand, gene_id))

        merged_exons = merge_intervals(exons.get(gene_id, []))
        if not merged_exons:
            continue

        cursor = gene_start
        for exon_start, exon_end in merged_exons:
            if cursor < exon_start:
                introns_by_chr[chrom].append((cursor, exon_start - 1, strand, gene_id))
            cursor = max(cursor, exon_end + 1)

        if cursor <= gene_end:
            introns_by_chr[chrom].append((cursor, gene_end, strand, gene_id))

    return introns_by_chr, gene_spans_by_chr

### COVERAGE-BASED UNIQUENESS CALCULATION #####################################


def build_segments(events: list[tuple[int, int]]) -> tuple[list[int], list[int]]:
    """Build stepwise coverage segments from sweep-line events."""

    points: list[int] = []
    coverage: list[int] = []
    active = 0
    index = 0

    while index < len(events):
        position = events[index][0]
        points.append(position)
        while index < len(events) and events[index][0] == position:
            active += events[index][1]
            index += 1
        coverage.append(active)

    return points, coverage


def bases_covered_by_exactly_one_gene(
    start: int,
    end: int,
    points: list[int],
    coverage: list[int],
) -> int:
    """Count bases with coverage == 1 across [start, end]."""

    index = bisect_left(points, start)
    if index == len(points) or points[index] > start:
        index -= 1
    index = max(index, 0)

    total = 0
    while index < len(points):
        segment_start = max(start, points[index])
        segment_end = end
        if index + 1 < len(points):
            segment_end = min(end, points[index + 1] - 1)

        if segment_start <= segment_end and coverage[index] == 1:
            total += segment_end - segment_start + 1

        if index + 1 >= len(points) or points[index + 1] > end:
            break
        index += 1

    return total


def summarize_intron_uniqueness(
    introns_by_chr: dict[str, list[tuple[int, int, str, str]]],
    gene_spans_by_chr: dict[str, list[tuple[int, int, str, str]]],
) -> dict[str, float]:
    """Summarize how often annotated introns are unique to one gene."""
    total_intronic_intervals = 0
    total_intronic_bases = 0
    unique_same_strand_bases = 0
    unique_any_strand_bases = 0
    fully_unique_same_strand_intervals = 0
    fully_unique_any_strand_intervals = 0

    for chrom, introns in introns_by_chr.items():
        gene_spans = gene_spans_by_chr[chrom]

        events_any: list[tuple[int, int]] = []
        events_plus: list[tuple[int, int]] = []
        events_minus: list[tuple[int, int]] = []

        for gene_start, gene_end, strand, _gene_id in gene_spans:
            events_any.append((gene_start, 1))
            events_any.append((gene_end + 1, -1))
            if strand == "+":
                events_plus.append((gene_start, 1))
                events_plus.append((gene_end + 1, -1))
            else:
                events_minus.append((gene_start, 1))
                events_minus.append((gene_end + 1, -1))

        events_any.sort()
        events_plus.sort()
        events_minus.sort()

        points_any, cov_any = build_segments(events_any)
        points_plus, cov_plus = build_segments(events_plus)
        points_minus, cov_minus = build_segments(events_minus)

        for intron_start, intron_end, strand, _gene_id in introns:
            if intron_start > intron_end:
                continue

            intron_length = intron_end - intron_start + 1
            total_intronic_intervals += 1
            total_intronic_bases += intron_length

            if strand == "+":
                same_strand_unique = bases_covered_by_exactly_one_gene(
                    intron_start, intron_end, points_plus, cov_plus
                )
            else:
                same_strand_unique = bases_covered_by_exactly_one_gene(
                    intron_start, intron_end, points_minus, cov_minus
                )

            any_strand_unique = bases_covered_by_exactly_one_gene(
                intron_start, intron_end, points_any, cov_any
            )

            unique_same_strand_bases += same_strand_unique
            unique_any_strand_bases += any_strand_unique

            if same_strand_unique == intron_length:
                fully_unique_same_strand_intervals += 1
            if any_strand_unique == intron_length:
                fully_unique_any_strand_intervals += 1

    return {
        "intronic_intervals": total_intronic_intervals,
        "intronic_bases_total": total_intronic_bases,
        "bases_unique_same_strand_pct": round(
            100 * unique_same_strand_bases / total_intronic_bases, 2
        ),
        "bases_unique_any_strand_pct": round(
            100 * unique_any_strand_bases / total_intronic_bases, 2
        ),
        "intervals_fully_unique_same_strand_pct": round(
            100 * fully_unique_same_strand_intervals / total_intronic_intervals, 2
        ),
        "intervals_fully_unique_any_strand_pct": round(
            100 * fully_unique_any_strand_intervals / total_intronic_intervals, 2
        ),
    }

### OUTPUT WRITING ############################################################


def write_metrics_csv(metrics: dict[str, float], output_path: Path) -> None:
    """Write the final metrics as a simple two-column CSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        for metric, value in metrics.items():
            writer.writerow([metric, value])

### MAIN WORKFLOW #############################################################


def main() -> None:
    """Run the full intron-uniqueness calculation and report the result."""
    args = parse_args()

    genes, exons = load_gene_and_exon_annotations(args.gtf)
    introns_by_chr, gene_spans_by_chr = derive_introns(genes, exons)
    metrics = summarize_intron_uniqueness(introns_by_chr, gene_spans_by_chr)
    write_metrics_csv(metrics, args.output)

    print(f"Analyzed GTF: {args.gtf}")
    print(f"Wrote metrics: {args.output}")
    for metric, value in metrics.items():
        print(f"{metric},{value}")


if __name__ == "__main__":
    main()
