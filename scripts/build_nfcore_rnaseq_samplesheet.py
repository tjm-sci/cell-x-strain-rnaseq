#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
from datetime import datetime
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_METADATA = REPO_ROOT / "metadata" / "25-04-07-EXP_383_FANS_SAMPLE_METADATA.csv"
MERGED_FASTQ_MANIFEST = REPO_ROOT / "metadata" / "exp383_merged_fastq_manifest.csv"
OUTPUT_SAMPLESHEET = REPO_ROOT / "metadata" / "exp383_nfcore_rnaseq_samplesheet_template.csv"
MISSING_SAMPLES_CSV = REPO_ROOT / "metadata" / "MISSING.csv"
BATCH_SAMPLESHEET_DIR = REPO_ROOT / "metadata" / "batches"

POPULATIONS = ["NeuN", "SOX10", "SOX2", "PU1"]

OUTPUT_COLUMNS = [
    "sample",
    "fastq_1",
    "fastq_2",
    "strandedness",
    "mouse_n",
    "population",
    "inoculum",
    "dpi",
    "group_assignment",
    "inoculation_batch",
    "sample_mass_mg",
    "date_nuc_prep",
    "incubation_time_hrs",
]


def parse_group_assignment(group_assignment: str) -> tuple[str, str]:
    match = re.fullmatch(r"group_([A-Za-z0-9]+)_(\d+)", group_assignment)
    if not match:
        raise ValueError(f"Unexpected group_assignment value: {group_assignment}")
    inoculum, dpi = match.groups()
    return inoculum, dpi


def format_prep_date(raw_date: str) -> str:
    return datetime.strptime(raw_date, "%y%m%d").strftime("%Y-%m-%d")


def build_sample_id(mouse_n: str, inoculum: str, dpi: str, group_no: str, population: str) -> str:
    return f"{mouse_n}-{inoculum}-{dpi}-{group_no}-{population}"


def load_missing_sample_ids() -> set[str]:
    with MISSING_SAMPLES_CSV.open(newline="") as handle:
        return {row["sample"] for row in csv.DictReader(handle) if row.get("sample")}


def load_merged_fastq_paths() -> dict[str, dict[str, str]]:
    with MERGED_FASTQ_MANIFEST.open(newline="") as handle:
        manifest_rows = list(csv.DictReader(handle))

    fastqs: dict[str, dict[str, str]] = {}
    for row in manifest_rows:
        sample = row["sample"]
        if sample in fastqs:
            raise ValueError(f"Duplicate sample in merged FASTQ manifest: {sample}")
        fastq_1 = row["fastq_1"]
        fastq_2 = row["fastq_2"]
        if not fastq_1 or not fastq_2:
            raise ValueError(f"Missing merged FASTQ path(s) for sample: {sample}")
        fastqs[sample] = {"fastq_1": fastq_1, "fastq_2": fastq_2}

    return fastqs


def build_sample_rows() -> list[dict[str, str]]:
    with SOURCE_METADATA.open(newline="") as handle:
        source_rows = list(csv.DictReader(handle))

    if len(source_rows) != 72:
        raise ValueError(f"Expected 72 mouse rows, found {len(source_rows)}")

    merged_fastqs = load_merged_fastq_paths()
    sample_rows: list[dict[str, str]] = []
    seen_samples: set[str] = set()
    missing_sample_ids = load_missing_sample_ids()

    for source_row in source_rows:
        inoculum, dpi = parse_group_assignment(source_row["group_assignment"])
        prep_date = format_prep_date(source_row["date_nuc_prep"])

        for population in POPULATIONS:
            sample_id = build_sample_id(
                mouse_n=source_row["mouse_n"],
                inoculum=inoculum,
                dpi=dpi,
                group_no=source_row["group_no"],
                population=population,
            )
            if sample_id in missing_sample_ids:
                continue
            if sample_id in seen_samples:
                raise ValueError(f"Duplicate sample name generated: {sample_id}")
            seen_samples.add(sample_id)
            if sample_id not in merged_fastqs:
                raise ValueError(f"Sample missing from merged FASTQ manifest: {sample_id}")

            sample_rows.append(
                {
                    "sample": sample_id,
                    "fastq_1": merged_fastqs[sample_id]["fastq_1"],
                    "fastq_2": merged_fastqs[sample_id]["fastq_2"],
                    "strandedness": "auto",
                    "mouse_n": source_row["mouse_n"],
                    "population": population,
                    "inoculum": inoculum,
                    "dpi": dpi,
                    "group_assignment": source_row["group_assignment"],
                    "inoculation_batch": source_row["inoculation_batch"],
                    "sample_mass_mg": source_row["sample_mass_mg"],
                    "date_nuc_prep": prep_date,
                    "incubation_time_hrs": source_row["incubation_time_hrs"],
                }
            )

    expected_rows = 72 * 4 - len(missing_sample_ids)
    if len(sample_rows) != expected_rows:
        raise ValueError(f"Expected {expected_rows} RNA-seq rows, found {len(sample_rows)}")
    extra_manifest_samples = sorted(set(merged_fastqs) - seen_samples)
    if extra_manifest_samples:
        raise ValueError(
            f"Merged FASTQ manifest contains samples absent from metadata-derived sheet: {extra_manifest_samples[:10]}"
        )

    return sample_rows


def write_samplesheet(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def write_population_batches(rows: list[dict[str, str]]) -> None:
    BATCH_SAMPLESHEET_DIR.mkdir(parents=True, exist_ok=True)
    for population in POPULATIONS:
        population_rows = [row for row in rows if row["population"] == population]
        if not population_rows:
            raise ValueError(f"No rows found for population: {population}")
        output_path = BATCH_SAMPLESHEET_DIR / f"exp383_nfcore_rnaseq_{population}.csv"
        write_samplesheet(output_path, population_rows)
        print(f"Wrote {len(population_rows)} rows to {output_path}")


def main() -> None:
    rows = build_sample_rows()
    write_samplesheet(OUTPUT_SAMPLESHEET, rows)
    print(f"Wrote {len(rows)} rows to {OUTPUT_SAMPLESHEET}")
    write_population_batches(rows)


if __name__ == "__main__":
    main()
