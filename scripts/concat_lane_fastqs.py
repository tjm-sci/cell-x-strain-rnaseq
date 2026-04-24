#!/usr/bin/env python3
"""Concatenate lane-level FASTQs into one merged R1/R2 pair per sample.

This script exists because the raw data originally arrived lane-split, while
the storage strategy for the main nf-core run worked better with one merged
FASTQ pair per biological sample.
"""

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

### COMMAND-LINE ARGUMENTS ####################################################


def parse_args():
    """Parse the small set of inputs needed for one merge pass."""
    p = argparse.ArgumentParser(description='Concatenate lane-level gzipped FASTQs into merged per-sample FASTQs.')
    p.add_argument('--samplesheet', required=True, help='Lane-wise nf-core samplesheet CSV')
    p.add_argument('--output-dir', required=True, help='Directory for merged FASTQs')
    p.add_argument('--manifest', required=True, help='Output manifest CSV recording merged files')
    p.add_argument('--resume', action='store_true', help='Skip samples with existing non-empty merged outputs')
    return p.parse_args()

### HELPER FUNCTIONS ##########################################################


def lane_sort_key(row):
    """Sort lanes deterministically before concatenation."""
    lane = row.get('lane', '') or ''
    return (lane, row['fastq_1'], row['fastq_2'])


def ensure_exists(paths):
    """Fail early if any expected lane-level FASTQ is missing."""
    missing = [p for p in paths if not Path(p).exists()]
    if missing:
        for p in missing[:10]:
            print(f'MISSING\t{p}', file=sys.stderr)
        raise FileNotFoundError(f'{len(missing)} input FASTQ paths do not exist')


def concatenate(inputs, output_path):
    """Concatenate gzipped FASTQs byte-for-byte without decompressing."""
    tmp_path = output_path.with_suffix(output_path.suffix + '.tmp')
    with tmp_path.open('wb') as out_handle:
        subprocess.run(['cat', *inputs], check=True, stdout=out_handle)
    tmp_path.replace(output_path)

### MAIN WORKFLOW #############################################################


def main():
    """Group lanes by sample, concatenate them, and write a merge manifest."""
    args = parse_args()
    samplesheet = Path(args.samplesheet)
    output_dir = Path(args.output_dir)
    manifest_path = Path(args.manifest)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    grouped = defaultdict(list)
    with samplesheet.open(newline='') as handle:
        reader = csv.DictReader(handle)
        required = ['sample', 'fastq_1', 'fastq_2']
        missing_cols = [c for c in required if c not in reader.fieldnames]
        if missing_cols:
            raise ValueError(f'Samplesheet missing required columns: {missing_cols}')
        for row in reader:
            grouped[row['sample']].append(row)

    rows_out = []
    total = len(grouped)
    for idx, sample in enumerate(sorted(grouped), start=1):
        rows = sorted(grouped[sample], key=lane_sort_key)
        r1_inputs = [row['fastq_1'] for row in rows]
        r2_inputs = [row['fastq_2'] for row in rows]
        ensure_exists(r1_inputs + r2_inputs)

        out_r1 = output_dir / f'{sample}_R1.fastq.gz'
        out_r2 = output_dir / f'{sample}_R2.fastq.gz'

        if args.resume and out_r1.exists() and out_r2.exists() and out_r1.stat().st_size > 0 and out_r2.stat().st_size > 0:
            status = 'skipped_existing'
        else:
            print(f'[{idx}/{total}] merging {sample} ({len(rows)} lanes)', flush=True)
            concatenate(r1_inputs, out_r1)
            concatenate(r2_inputs, out_r2)
            status = 'merged'

        rows_out.append({
            'sample': sample,
            'fastq_1': str(out_r1),
            'fastq_2': str(out_r2),
            'lane_count': len(rows),
            'status': status,
        })

    with manifest_path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['sample', 'fastq_1', 'fastq_2', 'lane_count', 'status'])
        writer.writeheader()
        writer.writerows(rows_out)

    print(f'Wrote manifest: {manifest_path}', flush=True)


if __name__ == '__main__':
    main()
