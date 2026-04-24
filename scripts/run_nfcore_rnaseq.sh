#!/usr/bin/env bash
set -euo pipefail

### SCRIPT SETUP ##############################################################

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

export CONDARC="${repo_root}/config/condarc.nfcore.yml"
export NXF_HOME="${NXF_HOME:-/media/tmurphy/2TB_SSD/exp383/.nextflow}"
export NXF_CONDA_CACHEDIR="${NXF_CONDA_CACHEDIR:-/media/tmurphy/2TB_SSD/exp383/nf_conda_cache}"

pipeline_version="3.23.0"
params_file="${repo_root}/params/nfcore_rnaseq_fastp_salmon_rrna.yaml"
batch_sheet_dir="${repo_root}/metadata/batches"
results_root="/media/tmurphy/4TB_HDD/exp383/nfcore_rnaseq"
work_root="/media/tmurphy/2TB_SSD/exp383/nfcore_rnaseq_work"
logs_root="${repo_root}/logs/nextflow/batches"
shared_reference_root="/media/tmurphy/4TB_HDD/exp383/reference_cache"

populations=(NeuN SOX10 SOX2 PU1)
selected_populations=()
nextflow_args=()
cleanup_work="true"

### COMMAND-LINE ARGUMENTS ####################################################

# Keep the wrapper intentionally narrow:
# - optional per-population reruns
# - optional work-directory retention
# - everything else passes straight through to Nextflow
while (($#)); do
  case "$1" in
    --population)
      selected_populations+=("$2")
      shift 2
      ;;
    --population=*)
      selected_populations+=("${1#*=}")
      shift
      ;;
    --keep-work)
      cleanup_work="false"
      shift
      ;;
    *)
      nextflow_args+=("$1")
      shift
      ;;
    esac
done

### POPULATION SELECTION ######################################################

if ((${#selected_populations[@]} == 0)); then
  selected_populations=("${populations[@]}")
fi

ordered_populations=()
for population in "${populations[@]}"; do
  for requested in "${selected_populations[@]}"; do
    if [[ "${requested}" == "${population}" ]]; then
      ordered_populations+=("${population}")
    fi
  done
done

if ((${#ordered_populations[@]} == 0)); then
  echo "No valid populations selected. Choose from: ${populations[*]}" >&2
  exit 1
fi

### SHARED DIRECTORY SETUP ####################################################

mkdir -p "${logs_root}" "${results_root}" "${work_root}" "${shared_reference_root}" "${NXF_HOME}" "${NXF_CONDA_CACHEDIR}"

echo "Rebuilding merged master samplesheet and population batch sheets..."
python3 "${repo_root}/scripts/build_nfcore_rnaseq_samplesheet.py"
echo

total_runs="${#ordered_populations[@]}"
run_index=0

### BATCH LOOP ################################################################

# Each sorted population is run separately so that the working SSD can be
# recycled between batches.
for population in "${ordered_populations[@]}"; do
  run_index="$((run_index + 1))"
  case "${population}" in
    NeuN) batch_prefix="01" ;;
    SOX10) batch_prefix="02" ;;
    SOX2) batch_prefix="03" ;;
    PU1) batch_prefix="04" ;;
    *)
      echo "Unknown population tag mapping for: ${population}" >&2
      exit 1
      ;;
  esac
  batch_tag="${batch_prefix}_${population}"
  input_sheet="${batch_sheet_dir}/exp383_nfcore_rnaseq_${population}.csv"
  outdir="${results_root}/${batch_tag}"
  workdir="${work_root}/${batch_tag}"
  batch_log_dir="${logs_root}/${batch_tag}"
  console_log="${batch_log_dir}/console.log"
  nextflow_log="${batch_log_dir}/nextflow.log"
  report_file="${batch_log_dir}/execution_report.html"
  trace_file="${batch_log_dir}/execution_trace.tsv"
  timeline_file="${batch_log_dir}/execution_timeline.html"
  dag_file="${batch_log_dir}/pipeline_dag.svg"
  salmon_index_dir="${shared_reference_root}/salmon"

  mkdir -p "${batch_log_dir}" "${outdir}" "${workdir}"
  export NXF_WORK="${workdir}"

  if [[ ! -f "${input_sheet}" ]]; then
    echo "Missing batch samplesheet: ${input_sheet}" >&2
    exit 1
  fi

  sample_count="$(( $(wc -l < "${input_sheet}") - 1 ))"
  rm -f "${console_log}" "${nextflow_log}" "${report_file}" "${trace_file}" "${timeline_file}" "${dag_file}"

  # Build the Nextflow command as an array so additional user-supplied args can
  # be appended safely without shell quoting problems.
  nf_cmd=(
    mamba run -n rnaseq
    nextflow
    -log "${nextflow_log}"
    run nf-core/rnaseq
    -r "${pipeline_version}"
    -profile conda
    -c "${repo_root}/nextflow.config"
    -params-file "${params_file}"
    --input "${input_sheet}"
    --outdir "${outdir}"
    -ansi-log false
    -with-report "${report_file}"
    -with-trace "${trace_file}"
    -with-timeline "${timeline_file}"
    -with-dag "${dag_file}"
  )

  if [[ -d "${salmon_index_dir}" ]]; then
    nf_cmd+=(--salmon_index "${salmon_index_dir}")
  fi

  if ((${#nextflow_args[@]} > 0)); then
    nf_cmd+=("${nextflow_args[@]}")
  fi

  ### BATCH EXECUTION #########################################################

  printf '\n============================================================\n'
  printf 'Batch %d/%d\n' "${run_index}" "${total_runs}"
  printf 'Population: %s\n' "${population}"
  printf 'Samplesheet: %s\n' "${input_sheet}"
  printf 'Samples: %s\n' "${sample_count}"
  printf 'Results dir: %s\n' "${outdir}"
  printf 'Work dir: %s\n' "${workdir}"
  printf 'Console log: %s\n' "${console_log}"
  printf 'Nextflow log: %s\n' "${nextflow_log}"
  printf 'NXF_HOME: %s\n' "${NXF_HOME}"
  printf 'NXF_CONDA_CACHEDIR: %s\n' "${NXF_CONDA_CACHEDIR}"
  if [[ -d "${salmon_index_dir}" ]]; then
    printf 'Reusing Salmon index: %s\n' "${salmon_index_dir}"
  else
    printf 'Reusing Salmon index: no\n'
  fi
  printf 'Command:\n'
  printf '  %q' "${nf_cmd[@]}"
  printf '\n============================================================\n\n'

  if ! "${nf_cmd[@]}" 2>&1 | tee "${console_log}"; then
    printf '\nBatch %d/%d failed for population %s.\n' "${run_index}" "${total_runs}" "${population}" >&2
    printf 'Work directory retained for inspection/resume: %s\n' "${workdir}" >&2
    exit 1
  fi

  if [[ -d "${outdir}/genome/index/salmon" ]]; then
    ln -sfn "${outdir}/genome/index/salmon" "${salmon_index_dir}"
    printf 'Registered reusable Salmon index: %s -> %s\n' "${salmon_index_dir}" "${outdir}/genome/index/salmon"
  fi

  if [[ "${cleanup_work}" == "true" ]]; then
    printf 'Deleting completed work directory to free SSD: %s\n' "${workdir}"
    rm -rf "${workdir}"
  else
    printf 'Keeping completed work directory: %s\n' "${workdir}"
  fi
done

### FINAL STATUS ##############################################################

printf '\nAll requested batches completed successfully.\n'
