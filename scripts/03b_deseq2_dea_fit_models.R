#!/usr/bin/env Rscript

# Fit population-wise DESeq2 models for EXP383.
#
# Inputs:
# - the label-corrected output object from 01e_sample_label_correction.R
# - the recommended population-specific design formulae from 02b
#
# Outputs:
# - one fitted DESeq2 object per population
# - one model-aware VST object per population
# - compact model QC tables and plots
#
# Deliberate boundary:
# This script fits models only. Contrasts, volcano plots, MA plots, and DEG
# tables are handled by 03c_deseq2_dea_export_results.R.

suppressPackageStartupMessages({
  required_packages <- c("DESeq2", "SummarizedExperiment", "ggplot2", "pheatmap", "here")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      ". Install them in the rnaseq environment before running this script.",
      call. = FALSE
    )
  }
})

suppressMessages(here::i_am("scripts/03b_deseq2_dea_fit_models.R"))
source(here::here("scripts", "path_helpers.R"))
source(here::here("scripts", "plot_style.R"))
source(here::here("scripts", "03a_deseq2_dea_utilities.R"))

# ----------------
# Fixed settings
# ----------------
# Defaults run all populations using the current corrected input and the
# population-specific design formulae chosen during covariate screening.
SETTINGS <- list(
  population = "all",
  input_rds = "results/dea/01e_sample_label_correction/exp383_salmon_gene_input_label_corrected.rds",
  screening_root = "results/dea/02_covariate_screening",
  output_root = "results/dea/03_deseq2_dea",
  min_count = "10",
  min_samples = "auto"
)

args <- parse_key_value_args(SETTINGS)

input_rds <- resolve_project_path(args$input_rds, must_work = TRUE)
screening_root <- resolve_project_path(args$screening_root, must_work = TRUE)
output_root <- resolve_project_path(args$output_root)

min_count <- as.integer(args$min_count)
if (is.na(min_count) || min_count < 0) {
  stop("--min_count must be a non-negative integer.", call. = FALSE)
}

min_samples <- if (identical(args$min_samples, "auto")) {
  "auto"
} else {
  as.integer(args$min_samples)
}
if (!identical(min_samples, "auto") && (is.na(min_samples) || min_samples < 1)) {
  stop("--min_samples must be 'auto' or a positive integer.", call. = FALSE)
}

# ----------------
# Load inputs
# ----------------
analysis_input <- load_dea_input(input_rds)
metadata <- prepare_coldata_for_deseq(analysis_input$sample_metadata)

available_populations <- order_known_levels(metadata$population, exp383_population_levels)
populations_to_run <- parse_population_selection(args$population, available_populations)

ensure_dir(output_root)

# -----------------------------
# Fit one model per population
# -----------------------------
run_summaries <- list()

for (population in populations_to_run) {
  message("Fitting DESeq2 model for ", population)

  # Carry forward the exact formula selected during covariate screening.
  design_formula <- read_population_design_formula(
    screening_root = screening_root,
    population = population
  )

  # Re-apply the same gene filter used for covariate screening so model input is
  # explicit and reproducible from this script alone.
  prepared <- prepare_population_for_deseq(
    analysis_input = analysis_input,
    population = population,
    design_formula = design_formula,
    min_count = min_count,
    min_samples = min_samples
  )

  # Stop early if the chosen formula is not estimable for this population.
  rank_check <- check_design_matrix_full_rank(
    design_formula = design_formula,
    coldata = prepared$coldata
  )
  if (!rank_check$is_full_rank) {
    stop(
      sprintf(
        "%s design matrix is not full rank: rank %s for %s columns.",
        population,
        rank_check$rank,
        rank_check$n_columns
      ),
      call. = FALSE
    )
  }

  population_root <- ensure_dir(file.path(output_root, population))
  model_root <- ensure_dir(file.path(population_root, "model_fit"))
  plot_dir <- ensure_dir(file.path(model_root, "plots"))

  # DESeqDataSetFromTximport uses the tximport abundance/length information to
  # preserve transcript-length correction while fitting at gene level.
  dds <- DESeq2::DESeqDataSetFromTximport(
    txi = prepared$txi,
    colData = prepared$coldata,
    design = stats::as.formula(design_formula)
  )

  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  vsd <- DESeq2::vst(dds, blind = FALSE)

  saveRDS(dds, file.path(model_root, "dds.rds"))
  saveRDS(vsd, file.path(model_root, "vst.rds"))

  # ----------------
  # Model metadata
  # ----------------
  writeLines(design_formula, file.path(model_root, "model_formula.txt"))
  writeLines(DESeq2::resultsNames(dds), file.path(model_root, "results_names.txt"))
  write_tsv(prepared$filter_summary, file.path(model_root, "filter_summary.tsv"))

  design_matrix <- as.data.frame(rank_check$design_matrix, check.names = FALSE)
  design_matrix$sample <- rownames(rank_check$design_matrix)
  design_matrix <- design_matrix[, c("sample", setdiff(colnames(design_matrix), "sample")), drop = FALSE]
  write_tsv(design_matrix, file.path(model_root, "design_matrix.tsv"))

  sample_qc <- as.data.frame(SummarizedExperiment::colData(dds))
  sample_qc$sample <- rownames(sample_qc)
  sample_qc$size_factor <- get_sample_scaling_factors(dds)
  sample_qc <- sample_qc[, c("sample", setdiff(colnames(sample_qc), "sample")), drop = FALSE]
  write_tsv(sample_qc, file.path(model_root, "sample_model_qc.tsv"))

  model_summary <- data.frame(
    population = population,
    n_samples = ncol(dds),
    n_genes = nrow(dds),
    design_formula = design_formula,
    design_rank = rank_check$rank,
    design_columns = rank_check$n_columns,
    min_count = min_count,
    min_samples = prepared$filter_summary$min_samples,
    stringsAsFactors = FALSE
  )
  write_tsv(model_summary, file.path(model_root, "model_fit_summary.tsv"))
  run_summaries[[population]] <- model_summary

  # ----------------
  # Model QC plots
  # ----------------
  pca_table <- compute_pca_table(vsd, SummarizedExperiment::colData(dds))
  write_tsv(pca_table, file.path(model_root, "vst_pca_scores.tsv"))

  save_model_pca_plot(
    pca_table = pca_table,
    color_var = "group_assignment",
    title_text = paste0(population, ": fitted-model VST PCA"),
    output_file = file.path(plot_dir, "vst_pca_by_group_assignment.png")
  )

  save_sample_distance_heatmap(
    vsd = vsd,
    coldata = SummarizedExperiment::colData(dds),
    output_file = file.path(plot_dir, "sample_distance_heatmap.png"),
    population = population
  )

  save_size_factor_plot(
    size_factor_table = sample_qc,
    output_file = file.path(plot_dir, "size_factors.png"),
    population = population
  )

  save_dispersion_plot(
    dds = dds,
    output_file = file.path(plot_dir, "dispersion_estimates.png"),
    population = population
  )
}

combined_summary <- do.call(rbind, run_summaries)
rownames(combined_summary) <- NULL
write_tsv(combined_summary, file.path(output_root, "model_fit_summary.tsv"))

message("Finished DESeq2 model fitting.")
