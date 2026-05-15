#!/usr/bin/env Rscript

# This script implements two linked tasks for one sorted nuclei population per script run:
#
# 1. Formal technical covariate selection
#    We screen technical variables using:
#    - pairwise collinearity checks
#    - blind (no covariates included) PCA on the VST-tranformed count matrix
#    - weighted covariate-to-PC association metrics
#    - multivariable variancePartition
#    The goal is to return an explicit, machine-readable set of technical
#    covariates to carry forward into DEA.
#
# 2. DESeq2 design formula testing
#    We then fit a small number of interpretable DESeq2 formulas, including a
#    recommended formula built from the selected covariates, and write the usual
#    PCA / distance / size-factor / dispersion diagnostics for each.
#    The final formula will be used for DEA.
#
# The script requires the output object from 01_build_salmon_gene_inputs.R
# and the EXP383 metadata structure.

suppressPackageStartupMessages({
  required_packages <- c("DESeq2", "ggplot2", "pheatmap", "variancePartition", "here")
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

suppressMessages(here::i_am("scripts/02b_covariate_screening_main_script.R"))
source(here::here("scripts", "path_helpers.R"))

### Fixed settings ###
# Set paths and thresholds for the covariate screening and design formula testing. These are kept the same across populations.
SETTINGS <- list(
  input_rds = "results/dea/01_build_salmon_gene_inputs/exp383_salmon_gene_input.rds",
  output_root = "results/dea/02_covariate_screening",
  min_count = 10,
  min_samples = "auto",
  pc_variance_threshold = 0.85,
  pc_padj_cutoff = 0.05,
  weighted_pc_cutoff = 0.025,
  varpart_q3_cutoff = 0.03,
  varpart_max_cutoff = 0.75
)

# Load the shared plotting helpers from this repo so that DEA diagnostics and
# global QC figures use the same population colours and base theme.
source(here::here("scripts", "plot_style.R"))
source(here::here("scripts", "02a_covariate_screening_utilities.R"))

# -----------------
# Main script logic
# -----------------
args <- parse_cli_args()
SETTINGS$input_rds <- resolve_project_path(SETTINGS$input_rds)
SETTINGS$output_root <- resolve_project_path(SETTINGS$output_root)
input_rds <- SETTINGS$input_rds
design_tsv <- resolve_project_path(args$design_tsv)
stop_if_missing(input_rds, "Prepared input RDS")
stop_if_missing(design_tsv, "Design TSV")

analysis_input <- readRDS(input_rds)
required_objects <- c("txi_gene", "sample_metadata")
missing_objects <- setdiff(required_objects, names(analysis_input))
if (length(missing_objects) > 0) {
  stop(
    "Input RDS does not look like the output of 01_build_salmon_gene_inputs.R. ",
    "Missing objects: ", paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

txi_gene <- analysis_input$txi_gene
coldata_all <- analysis_input$sample_metadata
coldata_all <- prepare_coldata_for_modeling(coldata_all)
coldata_all <- ensure_standardized_covariates(coldata_all)
if ("population" %in% colnames(coldata_all)) {
  coldata_all$population <- exp383_population_factor(coldata_all$population)
}

if (!args$population %in% as.character(unique(coldata_all$population))) {
  stop(
    sprintf("Population '%s' not found in prepared metadata.", args$population),
    call. = FALSE
  )
}

keep_samples <- coldata_all$population == args$population
coldata_pop <- coldata_all[keep_samples, , drop = FALSE]
rownames(coldata_pop) <- coldata_pop$sample

txi_pop <- subset_txi(
  txi = txi_gene,
  keep_rows = rep(TRUE, nrow(txi_gene$counts)),
  keep_cols = keep_samples
)

if (identical(SETTINGS$min_samples, "auto")) {
  group_sizes <- table(coldata_pop$group_assignment)
  min_samples <- as.integer(min(group_sizes))
} else {
  min_samples <- SETTINGS$min_samples
}

keep_genes <- compute_shared_filter(
  count_matrix = txi_pop$counts,
  grouping = coldata_pop$group_assignment,
  min_count = SETTINGS$min_count,
  min_samples = min_samples
)

txi_pop_filtered <- subset_txi(
  txi = txi_pop,
  keep_rows = keep_genes,
  keep_cols = rep(TRUE, ncol(txi_pop$counts))
)

shared_root <- file.path(SETTINGS$output_root, args$population, "shared")
dir.create(shared_root, recursive = TRUE, showWarnings = FALSE)

filter_summary <- data.frame(
  population = args$population,
  n_samples = nrow(coldata_pop),
  n_genes_before_filter = nrow(txi_pop$counts),
  n_genes_after_filter = nrow(txi_pop_filtered$counts),
  min_count = SETTINGS$min_count,
  min_samples = min_samples,
  filter_group_var = "group_assignment",
  stringsAsFactors = FALSE
)

write_tsv(filter_summary, file.path(shared_root, "filter_summary.tsv"))

# ------------------------------
# Shared formal selection stage
# ------------------------------
message("Running shared covariate screening for population: ", args$population)

screen_root <- file.path(shared_root, "covariate_screening")
screen_plot_dir <- file.path(screen_root, "plots")
dir.create(screen_plot_dir, recursive = TRUE, showWarnings = FALSE)

registry <- build_covariate_registry(coldata_pop)
pairwise_associations <- build_pairwise_covariate_associations(coldata_pop, registry)
numeric_collinearity_summary <- build_numeric_collinearity_summary(pairwise_associations, registry)
formal_registry <- finalize_covariate_registry(registry, numeric_collinearity_summary, n_samples = nrow(coldata_pop))

dds_blind <- DESeq2::DESeqDataSetFromTximport(
  txi = txi_pop_filtered,
  colData = coldata_pop,
  design = stats::as.formula("~ 1")
)
dds_blind <- DESeq2::estimateSizeFactors(dds_blind)
vsd_blind <- DESeq2::vst(dds_blind, blind = TRUE)

blind_pca_table <- compute_pca_table(vsd_blind, coldata_pop)
percent_var <- attr(blind_pca_table, "percent_var")
selected_pcs <- select_pcs_by_cumulative_variance(
  percent_var,
  cumulative_threshold = SETTINGS$pc_variance_threshold
)

pc_variance_table <- data.frame(
  pc = paste0("PC", seq_along(percent_var)),
  pc_index = seq_along(percent_var),
  variance_percent = round(percent_var, 4),
  cumulative_variance_percent = round(cumsum(percent_var), 4),
  selected_for_screen = paste0("PC", seq_along(percent_var)) %in% selected_pcs,
  stringsAsFactors = FALSE
)

pc_covariate_associations <- build_pc_covariate_associations(
  pc_table = blind_pca_table,
  registry = formal_registry,
  selected_pcs = selected_pcs
)
pc_covariate_summary <- summarize_pc_covariate_associations(
  pc_association_table = pc_covariate_associations,
  registry = formal_registry,
  pc_padj_cutoff = SETTINGS$pc_padj_cutoff
)

variance_partition_results <- run_multivariate_variance_partition(
  vsd = vsd_blind,
  registry = formal_registry,
  coldata = coldata_pop
)
variance_partition_summary <- summarize_variance_partition(
  varpart_object = variance_partition_results$varpart,
  registry = formal_registry
)

selection_metrics <- apply_formal_selection_rule(
  registry = formal_registry,
  pc_summary = pc_covariate_summary,
  varpart_summary = variance_partition_summary,
  settings = SETTINGS
)
selected_covariates <- selection_metrics[selection_metrics$selected, , drop = FALSE]
recommended_formula <- build_recommended_design_formula(selection_metrics)

# Only the essential outputs are kept here: the final decision table, the
# selected covariates, the recommended design, and the design list that will be
# run in the second half of the script.
write_tsv(selection_metrics, file.path(screen_root, "covariate_selection_metrics.tsv"))
write_tsv(selected_covariates, file.path(screen_root, "selected_covariates.tsv"))
write_tsv(
  data.frame(
    pc_variance_threshold = SETTINGS$pc_variance_threshold,
    pc_padj_cutoff = SETTINGS$pc_padj_cutoff,
    weighted_pc_cutoff = SETTINGS$weighted_pc_cutoff,
    varpart_q3_cutoff = SETTINGS$varpart_q3_cutoff,
    varpart_max_cutoff = SETTINGS$varpart_max_cutoff,
    stringsAsFactors = FALSE
  ),
  file.path(screen_root, "covariate_selection_thresholds.tsv")
)
writeLines(recommended_formula, con = file.path(screen_root, "recommended_design_formula.txt"))

save_pca_scree_plot(
  pc_variance_table = pc_variance_table,
  n_selected_pcs = length(selected_pcs),
  variance_threshold = SETTINGS$pc_variance_threshold,
  output_file = file.path(screen_plot_dir, "blind_pca_scree.png")
)

save_numeric_correlation_heatmap(
  coldata = coldata_pop,
  registry = formal_registry,
  output_file = file.path(screen_plot_dir, "retained_numeric_covariate_correlations.png")
)

plot_selection_metrics <- selection_metrics[selection_metrics$is_candidate, , drop = FALSE]
save_metric_barplot(
  plot_table = plot_selection_metrics[!is.na(plot_selection_metrics$max_weighted_pc_score), , drop = FALSE],
  metric_column = "max_weighted_pc_score",
  title_text = "Maximum weighted PC association score",
  y_label = "Max weighted PC score",
  output_file = file.path(screen_plot_dir, "weighted_pc_scores.png")
)
save_metric_barplot(
  plot_table = plot_selection_metrics,
  metric_column = "varpart_q3_percent",
  title_text = "variancePartition Q3 by covariate",
  y_label = "Q3 variance explained (%)",
  output_file = file.path(screen_plot_dir, "variance_partition_q3.png")
)
save_metric_barplot(
  plot_table = plot_selection_metrics,
  metric_column = "varpart_max_percent",
  title_text = "variancePartition max by covariate",
  y_label = "Max variance explained (%)",
  output_file = file.path(screen_plot_dir, "variance_partition_max.png")
)
save_variance_partition_plot(
  varpart_object = variance_partition_results$varpart,
  output_file = file.path(screen_plot_dir, "variance_partition_multivariate.png")
)

blind_pca_plot_table <- blind_pca_table
attr(blind_pca_plot_table, "percent_var") <- percent_var[1:2]
shared_plot_covariates <- unique(c(
  "group_assignment",
  selection_metrics$covariate[selection_metrics$selected]
))
for (var_name in shared_plot_covariates) {
  save_pca_plot(
    pca_table = blind_pca_plot_table,
    color_var = var_name,
    title_text = sprintf("Blind PCA: %s", var_name),
    output_file = file.path(screen_plot_dir, sprintf("blind_pca_by_%s.png", var_name))
  )
}

# ---------------------------------------
# Design-specific DESeq2 QC and behaviour
# ---------------------------------------
design_table <- load_design_table(design_tsv)
design_table <- append_recommended_design(design_table, recommended_formula)
write_tsv(design_table, file.path(shared_root, "designs_run.tsv"))

for (i in seq_len(nrow(design_table))) {
  design_id <- design_table$design_id[[i]]
  design_formula <- design_table$design_formula[[i]]
  design_formula_obj <- stats::as.formula(design_formula)

  design_dir <- file.path(SETTINGS$output_root, args$population, design_id)
  plot_dir <- file.path(design_dir, "plots")
  table_dir <- file.path(design_dir, "tables")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  message("Testing design: ", design_id, " -> ", design_formula)

  writeLines(
    c(
      paste0("population: ", args$population),
      paste0("design_id: ", design_id),
      paste0("design_formula: ", design_formula),
      paste0("min_count: ", SETTINGS$min_count),
      paste0("min_samples: ", min_samples)
    ),
    con = file.path(design_dir, "design_metadata.txt")
  )

  dds <- DESeq2::DESeqDataSetFromTximport(
    txi = txi_pop_filtered,
    colData = coldata_pop,
    design = design_formula_obj
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  vsd <- DESeq2::vst(dds, blind = FALSE)

  sample_scaling_factors <- get_sample_scaling_factors(dds)
  size_factor_table <- data.frame(
    sample = colnames(dds),
    size_factor = sample_scaling_factors,
    group_assignment = coldata_pop$group_assignment,
    stringsAsFactors = FALSE
  )

  pca_table <- compute_pca_table(vsd, coldata_pop)
  attr(pca_table, "percent_var") <- attr(pca_table, "percent_var")[1:2]

  design_matrix <- model.matrix(design_formula_obj, data = as.data.frame(SummarizedExperiment::colData(dds)))
  write_tsv(as.data.frame(design_matrix), file.path(table_dir, "design_matrix.tsv"), row_names = TRUE)
  write_tsv(size_factor_table, file.path(table_dir, "size_factors.tsv"))
  writeLines(DESeq2::resultsNames(dds), con = file.path(table_dir, "results_names.txt"))

  save_sample_distance_heatmap(
    vsd = vsd,
    annotation_df = as.data.frame(coldata_pop[, intersect(
      c("group_assignment", "inoculation_batch"),
      colnames(coldata_pop)
    ), drop = FALSE]),
    output_file = file.path(plot_dir, "sample_distance_heatmap.png")
  )
  save_dispersion_plot(dds, file.path(plot_dir, "dispersion_plot.png"))
  save_size_factor_plot(size_factor_table, file.path(plot_dir, "size_factors.png"))

  # The design-comparison stage only needs PCA plots for the biology of interest
  # and the terms actually present in the current design formula.
  variables_to_plot <- unique(c("group_assignment", all.vars(design_formula_obj)))

  for (var_name in variables_to_plot) {
    if (!var_name %in% colnames(pca_table)) {
      next
    }

    save_pca_plot(
      pca_table = pca_table,
      color_var = var_name,
      title_text = sprintf("PCA: %s | %s", design_id, var_name),
      output_file = file.path(plot_dir, sprintf("pca_by_%s.png", var_name))
    )
  }
}

message("Finished covariate screening workflow for population: ", args$population)
