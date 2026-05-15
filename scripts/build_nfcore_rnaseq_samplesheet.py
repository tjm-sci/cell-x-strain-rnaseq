#!/usr/bin/env python3
"""Build nf-core/rnaseq sample sheets from EXP383 metadata and FASTQ manifests.

The script writes one master sample sheet plus one per sorted nuclei population.
It is called by scripts/00_run_nfcore_rnaseq.sh before launching Nextflow.
"""

from __future__ import annotations

import csv
import re
from datetime import datetime
from pathlib import Path


### Paths and fixed settings #################################################

REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_METADATA = REPO_ROOT / "metadata" / "25-04-07-EXP_383_FANS_SAMPLE_METADATA.csv"
MERGED_FASTQ_MANIFEST = REPO_ROOT / "metadata" / "exp383_merged_fastq_manifest.csv"
OUTPUT_SAMPLESHEET = REPO_ROOT / "metadata" / "exp383_nfcore_rnaseq_samplesheet_template.csv"
MISSING_SAMPLES_CSV = REPO_ROOT / "metadata" / "MISSING.csv"
BATCH_SAMPLESHEET_DIR = REPO_ROOT / "metadata" / "batches"
MOUSE_EXCLUSION_FLAGS = (
    REPO_ROOT.parent
    / "exp383_mouse_metadata"
    / "output"
    / "exp383_mouse_exclusion_flags.csv"
)

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


### Parsing and formatting helpers ###########################################

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


### Input readers ############################################################

def load_missing_sample_ids() -> set[str]:
    """Return sample IDs listed as intentionally absent from sequencing input."""
    with MISSING_SAMPLES_CSV.open(newline="") as handle:
        return {row["sample"] for row in csv.DictReader(handle) if row.get("sample")}


def load_excluded_mouse_ids() -> set[str]:
    """Load project-level mouse exclusions that apply to RNA-seq sample sheets."""
    if not MOUSE_EXCLUSION_FLAGS.is_file():
        raise FileNotFoundError(
            "Project-level mouse exclusion flags not found: "
            f"{MOUSE_EXCLUSION_FLAGS}. Run exp383_mouse_metadata/scripts/05_check_non_experimental_culls.R first."
        )

    with MOUSE_EXCLUSION_FLAGS.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    required_columns = {"animal_id", "exclude_from_downstream", "exclusion_applies_to"}
    missing_columns = required_columns.difference(rows[0].keys()) if rows else required_columns
    if missing_columns:
        raise ValueError(f"Mouse exclusion flags missing columns: {sorted(missing_columns)}")

    excluded_mouse_ids: set[str] = set()
    for row in rows:
        exclude = row["exclude_from_downstream"].strip().lower() in {"true", "1"}
        applies = row["exclusion_applies_to"].strip()
        if exclude and applies in {"all_exp383_repos", REPO_ROOT.name}:
            excluded_mouse_ids.add(row["animal_id"].strip())

    return excluded_mouse_ids


def load_merged_fastq_paths() -> dict[str, dict[str, str]]:
    """Return merged FASTQ paths keyed by generated EXP383 sample ID."""
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


### Sample row construction ##################################################

def build_sample_rows() -> list[dict[str, str]]:
    """Build nf-core/rnaseq rows from mouse metadata crossed with populations."""
    with SOURCE_METADATA.open(newline="") as handle:
        source_rows = list(csv.DictReader(handle))

    if len(source_rows) != 72:
        raise ValueError(f"Expected 72 mouse rows, found {len(source_rows)}")

    # Remove mice culled for non-experimental reasons before expanding each
    # mouse into one row per sorted nuclei population.
    excluded_mouse_ids = load_excluded_mouse_ids()
    source_row_count = len(source_rows)
    source_rows = [
        row for row in source_rows
        if row["mouse_n"].strip() not in excluded_mouse_ids
    ]
    excluded_source_rows = source_row_count - len(source_rows)
    if excluded_source_rows:
        print(
            "Excluded RNA-seq source mice for non-experimental cull reasons: "
            f"{excluded_source_rows}"
        )

    merged_fastqs = load_merged_fastq_paths()
    sample_rows: list[dict[str, str]] = []
    seen_samples: set[str] = set()
    missing_sample_ids = load_missing_sample_ids()
    expected_rows = 0

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
            expected_rows += 1
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

    if len(sample_rows) != expected_rows:
        raise ValueError(f"Expected {expected_rows} RNA-seq rows, found {len(sample_rows)}")
    extra_manifest_samples = sorted(set(merged_fastqs) - seen_samples)
    unexpected_extra_manifest_samples = [
        sample for sample in extra_manifest_samples
        if sample.split("-", 1)[0] not in excluded_mouse_ids
    ]
    if unexpected_extra_manifest_samples:
        raise ValueError(
            "Merged FASTQ manifest contains samples absent from metadata-derived sheet: "
            f"{unexpected_extra_manifest_samples[:10]}"
        )

    return sample_rows


### Output writers ###########################################################

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


### Main workflow ############################################################

def main() -> None:
    rows = build_sample_rows()
    write_samplesheet(OUTPUT_SAMPLESHEET, rows)
    print(f"Wrote {len(rows)} rows to {OUTPUT_SAMPLESHEET}")
    write_population_batches(rows)


if __name__ == "__main__":
    main()
