#!/usr/bin/env python3
"""Summarize completed nf-core/rnaseq batches for the EXP383 project."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Iterable


BATCHES = [
    ("01_NeuN", "NeuN"),
    ("02_SOX10", "SOX10"),
    ("03_SOX2", "SOX2"),
    ("04_PU1", "PU1"),
]

### COMMAND-LINE ARGUMENTS ####################################################


def parse_args() -> argparse.Namespace:
    """Parse the batch-results root and summary output directory."""
    parser = argparse.ArgumentParser(
        description="Summarize completed nf-core/rnaseq batch outputs for EXP383."
    )
    parser.add_argument(
        "--results-root",
        default="/media/tmurphy/4TB_HDD/exp383/nfcore_rnaseq",
        help="Root directory containing per-population nf-core/rnaseq batch outputs.",
    )
    parser.add_argument(
        "--outdir",
        default="reports/nfcore_rnaseq_summary",
        help="Output directory for summary files.",
    )
    return parser.parse_args()

### SMALL FILE HELPERS ########################################################


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a small tabular MultiQC export into memory."""
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def is_sample_row(sample_name: str) -> bool:
    """Exclude per-read companion rows such as 'Sample Read 1'."""
    return not (sample_name.endswith(" Read 1") or sample_name.endswith(" Read 2"))


def safe_float(value: str | None) -> float | None:
    """Convert a numeric-looking string to float while preserving empties as None."""
    if value is None:
        return None
    value = value.strip()
    if not value:
        return None
    return float(value)


def mean(values: Iterable[float | None]) -> float | None:
    """Average a sequence of optional floats, skipping missing entries."""
    filtered = [value for value in values if value is not None]
    if not filtered:
        return None
    return sum(filtered) / len(filtered)


def fmt(value: float | None, digits: int = 2) -> str:
    """Format summary values consistently for the Markdown report."""
    if value is None:
        return "NA"
    return f"{value:.{digits}f}"

### BATCH SUMMARIES ###########################################################


def summarise_batch(results_root: Path, batch_tag: str, population: str) -> dict[str, object]:
    """Summarize one per-population nf-core result directory."""
    batch_dir = results_root / batch_tag
    data_dir = batch_dir / "multiqc" / "multiqc_report_data"
    general_stats = read_tsv(data_dir / "multiqc_general_stats.txt")
    sample_rows = [row for row in general_stats if is_sample_row(row["Sample"])]

    salmon_dirs = [
        path for path in (batch_dir / "salmon").iterdir() if path.is_dir() and path.name != "deseq2_qc"
    ]

    library_types = sorted(
        {
            (row.get("salmon-library_types") or "").strip()
            for row in sample_rows
            if (row.get("salmon-library_types") or "").strip()
        }
    )

    return {
        "batch_tag": batch_tag,
        "population": population,
        "sample_count": len(sample_rows),
        "salmon_dir_count": len(salmon_dirs),
        "library_types": library_types,
        "raw_read_pairs_m": mean(
            safe_float(row.get("fastqc_raw-total_sequences")) / 2
            if safe_float(row.get("fastqc_raw-total_sequences")) is not None
            else None
            for row in sample_rows
        ),
        "raw_gc_pct": mean(safe_float(row.get("fastqc_raw-percent_gc")) for row in sample_rows),
        "raw_dup_pct": mean(safe_float(row.get("fastqc_raw-percent_duplicates")) for row in sample_rows),
        "fastp_surviving_pct": mean(safe_float(row.get("fastp-pct_surviving")) for row in sample_rows),
        "fastp_adapter_pct": mean(safe_float(row.get("fastp-pct_adapter")) for row in sample_rows),
        "fastp_q30_pct": mean(
            safe_float(row.get("fastp-after_filtering_q30_rate")) for row in sample_rows
        ),
        "trimmed_read_length": mean(
            safe_float(row.get("fastqc_trimmed-avg_sequence_length")) for row in sample_rows
        ),
        "rrna_alignment_pct": mean(
            safe_float(row.get("bowtie2_rrna_removal-overall_alignment_rate")) for row in sample_rows
        ),
        "filtered_read_pairs_m": mean(
            safe_float(row.get("fastqc_filtered-total_sequences")) / 2
            if safe_float(row.get("fastqc_filtered-total_sequences")) is not None
            else None
            for row in sample_rows
        ),
        "filtered_dup_pct": mean(
            safe_float(row.get("fastqc_filtered-percent_duplicates")) for row in sample_rows
        ),
        "salmon_mapped_pct": mean(safe_float(row.get("salmon-percent_mapped")) for row in sample_rows),
        "salmon_mapped_fragments_m": mean(
            safe_float(row.get("salmon-num_mapped")) for row in sample_rows
        ),
        "compatible_fragment_ratio_pct": mean(
            safe_float(row.get("salmon-compatible_fragment_ratio")) for row in sample_rows
        ),
        "results_size_gb": batch_dir.stat().st_size,
    }


def summarise_overall(batch_summaries: list[dict[str, object]]) -> dict[str, object]:
    """Combine per-population summaries into one overall dataset summary."""
    overall = {
        "population": "ALL",
        "sample_count": sum(int(batch["sample_count"]) for batch in batch_summaries),
        "salmon_dir_count": sum(int(batch["salmon_dir_count"]) for batch in batch_summaries),
        "library_types": sorted(
            {
                library_type
                for batch in batch_summaries
                for library_type in batch["library_types"]  # type: ignore[index]
            }
        ),
    }
    for key in [
        "raw_read_pairs_m",
        "raw_gc_pct",
        "raw_dup_pct",
        "fastp_surviving_pct",
        "fastp_adapter_pct",
        "fastp_q30_pct",
        "trimmed_read_length",
        "rrna_alignment_pct",
        "filtered_read_pairs_m",
        "filtered_dup_pct",
        "salmon_mapped_pct",
        "salmon_mapped_fragments_m",
        "compatible_fragment_ratio_pct",
    ]:
        weighted_sum = 0.0
        total_samples = 0
        for batch in batch_summaries:
            value = batch[key]
            sample_count = int(batch["sample_count"])
            if value is None:
                continue
            weighted_sum += float(value) * sample_count
            total_samples += sample_count
        overall[key] = weighted_sum / total_samples if total_samples else None
    return overall

### OUTPUT WRITING ############################################################


def write_population_summary_tsv(out_path: Path, rows: list[dict[str, object]]) -> None:
    """Write the machine-readable per-population summary table."""
    fieldnames = [
        "population",
        "sample_count",
        "salmon_dir_count",
        "library_types",
        "raw_read_pairs_m",
        "raw_gc_pct",
        "raw_dup_pct",
        "fastp_surviving_pct",
        "fastp_adapter_pct",
        "fastp_q30_pct",
        "trimmed_read_length",
        "rrna_alignment_pct",
        "filtered_read_pairs_m",
        "filtered_dup_pct",
        "salmon_mapped_pct",
        "salmon_mapped_fragments_m",
        "compatible_fragment_ratio_pct",
    ]
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: (
                    ",".join(row["library_types"]) if key == "library_types" else row.get(key)
                )
                for key in fieldnames
            })


def markdown_table(rows: list[dict[str, object]]) -> str:
    """Render the main Markdown summary table."""
    header = (
        "| Population | Samples | Raw pairs/sample (M) | fastp survive % | "
        "Adapter % | Q30 % | rRNA align % | Filtered pairs/sample (M) | "
        "Salmon mapped % | Mapped fragments/sample (M) | Library type |\n"
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n"
    )
    body = []
    for row in rows:
        body.append(
            "| {population} | {sample_count} | {raw} | {survive} | {adapter} | {q30} | {rrna} | "
            "{filtered} | {salmon_pct} | {salmon_mapped} | {library_types} |".format(
                population=row["population"],
                sample_count=row["sample_count"],
                raw=fmt(row["raw_read_pairs_m"]),
                survive=fmt(row["fastp_surviving_pct"]),
                adapter=fmt(row["fastp_adapter_pct"]),
                q30=fmt(row["fastp_q30_pct"]),
                rrna=fmt(row["rrna_alignment_pct"]),
                filtered=fmt(row["filtered_read_pairs_m"]),
                salmon_pct=fmt(row["salmon_mapped_pct"]),
                salmon_mapped=fmt(row["salmon_mapped_fragments_m"]),
                library_types=", ".join(row["library_types"]) if row["library_types"] else "NA",
            )
        )
    return header + "\n".join(body)


def write_markdown_summary(out_path: Path, batch_summaries: list[dict[str, object]], overall: dict[str, object]) -> None:
    """Write the project-facing Markdown summary report."""
    lines = [
        "# EXP383 nf-core/rnaseq Summary",
        "",
        "## Overview",
        "",
        "- Pipeline completed successfully in 4 sequential population batches: `NeuN`, `SOX10`, `SOX2`, and `PU1`.",
        "- Processing route: `FQ_LINT -> FastQC(raw) -> fastp -> FastQC(trimmed) -> Bowtie2 rRNA removal -> FastQC(filtered) -> Salmon pseudoalignment -> tximport/tximeta -> DESeq2 QC -> MultiQC`.",
        "- Reference: `Ensembl release 115 / GRCm39`.",
        "- Total completed biological samples: `{}`.".format(overall["sample_count"]),
        "- Total Salmon quant directories produced: `{}`.".format(overall["salmon_dir_count"]),
        "- Salmon inferred library type(s): `{}`.".format(
            ", ".join(overall["library_types"]) if overall["library_types"] else "NA"
        ),
        "",
        "## Population Metrics",
        "",
        markdown_table(batch_summaries + [overall]),
        "",
        "## Step-Level Averages",
        "",
        "- Raw sequencing: mean `{}` million read pairs per sample, mean GC `{}`%, mean duplicate rate `{}`%.".format(
            fmt(overall["raw_read_pairs_m"]), fmt(overall["raw_gc_pct"]), fmt(overall["raw_dup_pct"])
        ),
        "- fastp trimming: mean survival `{}`%, mean adapter-trimmed fraction `{}`%, mean post-filter Q30 `{}`%, mean trimmed read length `{}` bp.".format(
            fmt(overall["fastp_surviving_pct"]),
            fmt(overall["fastp_adapter_pct"]),
            fmt(overall["fastp_q30_pct"]),
            fmt(overall["trimmed_read_length"]),
        ),
        "- Bowtie2 rRNA removal: mean alignment rate to the rRNA reference `{}`%, leaving `{}` million read pairs per sample after filtering.".format(
            fmt(overall["rrna_alignment_pct"]),
            fmt(overall["filtered_read_pairs_m"]),
        ),
        "- Salmon quantification: mean mapped fraction `{}`%, mean mapped fragments `{}` million per sample, mean compatible fragment ratio `{}`%.".format(
            fmt(overall["salmon_mapped_pct"]),
            fmt(overall["salmon_mapped_fragments_m"]),
            fmt(overall["compatible_fragment_ratio_pct"]),
        ),
        "",
        "## Notes",
        "",
        "- The per-population summary table is also available as TSV and JSON in this same directory.",
        "- `rRNA align %` is the Bowtie2 overall alignment rate against the supplied rRNA reference, used here as the depletion / residual rRNA estimate.",
        "- Read-count values are reported as estimated paired-end read pairs per sample by halving the combined read totals reported by MultiQC.",
    ]
    out_path.write_text("\n".join(lines) + "\n")

### MAIN WORKFLOW #############################################################


def main() -> None:
    """Summarize all completed batches and write TSV, JSON, and Markdown outputs."""
    args = parse_args()
    results_root = Path(args.results_root).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    batch_summaries = [summarise_batch(results_root, batch_tag, population) for batch_tag, population in BATCHES]
    overall = summarise_overall(batch_summaries)

    write_population_summary_tsv(outdir / "population_metrics.tsv", batch_summaries + [overall])
    (outdir / "population_metrics.json").write_text(
        json.dumps({"batches": batch_summaries, "overall": overall}, indent=2) + "\n"
    )
    write_markdown_summary(outdir / "pipeline_summary.md", batch_summaries, overall)

    print(f"Wrote {outdir / 'pipeline_summary.md'}")
    print(f"Wrote {outdir / 'population_metrics.tsv'}")
    print(f"Wrote {outdir / 'population_metrics.json'}")


if __name__ == "__main__":
    main()
