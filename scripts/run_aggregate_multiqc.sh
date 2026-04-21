#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
results_root="${1:-/media/tmurphy/4TB_HDD/exp383/nfcore_rnaseq}"
outdir="${2:-${repo_root}/reports/all_populations_multiqc}"
multiqc_config="${3:-${repo_root}/config/multiqc_aggregate_config.yml}"

find_multiqc_bin() {
  if command -v multiqc >/dev/null 2>&1; then
    command -v multiqc
    return 0
  fi

  local cached
  cached="$(find /media/tmurphy/2TB_SSD/exp383/nf_conda_cache -path '*/bin/multiqc' -type f 2>/dev/null | sort | head -n 1 || true)"
  if [[ -n "${cached}" ]]; then
    printf '%s\n' "${cached}"
    return 0
  fi

  return 1
}

multiqc_bin="$(find_multiqc_bin || true)"
if [[ -z "${multiqc_bin}" ]]; then
  echo "Could not locate a multiqc executable. Install MultiQC or keep the nf-core conda cache available." >&2
  exit 1
fi

if [[ ! -f "${multiqc_config}" ]]; then
  echo "Missing MultiQC config: ${multiqc_config}" >&2
  exit 1
fi

batch_dirs=(
  "${results_root}/01_NeuN"
  "${results_root}/02_SOX10"
  "${results_root}/03_SOX2"
  "${results_root}/04_PU1"
)

for batch_dir in "${batch_dirs[@]}"; do
  if [[ ! -d "${batch_dir}" ]]; then
    echo "Missing batch directory: ${batch_dir}" >&2
    exit 1
  fi
done

mkdir -p "${outdir}"

echo "Using MultiQC executable: ${multiqc_bin}"
echo "Using MultiQC config: ${multiqc_config}"
echo "Aggregating completed batch outputs into: ${outdir}"

"${multiqc_bin}" \
  --force \
  --ignore "*/multiqc/*" \
  --config "${multiqc_config}" \
  --outdir "${outdir}" \
  --filename "exp383_all_populations_multiqc_report.html" \
  --title "EXP383 RNA-seq: All Populations" \
  --comment "Aggregated across NeuN, SOX10, SOX2, and PU1 nf-core/rnaseq batches." \
  "${batch_dirs[@]}"

echo "Aggregate MultiQC report written to: ${outdir}/exp383_all_populations_multiqc_report.html"
