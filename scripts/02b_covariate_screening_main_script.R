#!/usr/bin/env Rscript

# This script runs formal technical covariate screening for one sorted nuclei
# population per script run.
#
# We screen technical variables using:
# - pairwise collinearity checks
# - blind (no covariates included) PCA on the VST-transformed count matrix
# - weighted covariate-to-PC association metrics
# - multivariable variancePartition
#
# The script requires the label-corrected output object from
# 01e_sample_label_correction.R and the EXP383 metadata structure.

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
# Set paths and thresholds for covariate screening. These are kept the same
# across populations.
SETTINGS <- list(
  input_rds = "results/dea/01e_sample_label_correction/exp383_salmon_gene_input_label_corrected.rds",
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
stop_if_missing(input_rds, "Label-corrected input RDS")

analysis_input <- readRDS(input_rds)
required_objects <- c("txi_gene", "sample_metadata")
missing_objects <- setdiff(required_objects, names(analysis_input))
if (length(missing_objects) > 0) {
  stop(
    "Input RDS does not look like the corrected output used for covariate screening. ",
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

population_root <- file.path(SETTINGS$output_root, args$population)
dir.create(population_root, recursive = TRUE, showWarnings = FALSE)

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

write_tsv(filter_summary, file.path(population_root, "filter_summary.tsv"))

# ----------------------------
# Formal covariate selection
# ----------------------------
message("Running covariate screening for population: ", args$population)

screen_root <- file.path(population_root, "covariate_screening")
screen_plot_dir <- file.path(screen_root, "plots")
biology_balance_plot_dir <- file.path(screen_plot_dir, "selected_covariates_by_biological_variables")
biology_balance_data_dir <- file.path(screen_root, "selected_covariates_by_biological_variables_data")
dir.create(screen_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(biology_balance_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(biology_balance_data_dir, recursive = TRUE, showWarnings = FALSE)

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

# These outputs check whether final selected technical covariates are balanced
# across the biological variables that define the planned DEA contrasts.
biology_association_table <- build_selected_covariate_biology_associations(
  coldata = coldata_pop,
  selection_metrics = selection_metrics
)
numeric_biology_plot_data <- build_selected_numeric_covariate_biology_data(
  coldata = coldata_pop,
  selection_metrics = selection_metrics
)
categorical_biology_plot_data <- build_selected_categorical_covariate_biology_data(
  coldata = coldata_pop,
  selection_metrics = selection_metrics
)

# Keep the output focused on the screen itself: the full decision table, selected
# covariates, thresholds, recommended design, and the plots that support the
# covariate-selection decision.
write_tsv(selection_metrics, file.path(screen_root, "covariate_selection_metrics.tsv"))
write_tsv(selected_covariates, file.path(screen_root, "selected_covariates.tsv"))
write_tsv(
  biology_association_table,
  file.path(biology_balance_data_dir, "selected_covariate_biological_associations.tsv")
)
write_tsv(
  numeric_biology_plot_data,
  file.path(biology_balance_data_dir, "selected_numeric_covariates_by_biological_variables.tsv")
)
write_tsv(
  categorical_biology_plot_data,
  file.path(biology_balance_data_dir, "selected_categorical_covariates_by_biological_variables.tsv")
)
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
  output_file = file.path(screen_plot_dir, "blind_pca_scree.png"),
  population = args$population
)

save_numeric_correlation_heatmap(
  coldata = coldata_pop,
  registry = formal_registry,
  output_file = file.path(screen_plot_dir, "retained_numeric_covariate_correlations.png"),
  population = args$population
)

plot_selection_metrics <- selection_metrics[selection_metrics$is_candidate, , drop = FALSE]
save_metric_barplot(
  plot_table = plot_selection_metrics[!is.na(plot_selection_metrics$max_weighted_pc_score), , drop = FALSE],
  metric_column = "max_weighted_pc_score",
  title_text = "Maximum weighted PC association score",
  y_label = "Max weighted PC score",
  output_file = file.path(screen_plot_dir, "weighted_pc_scores.png"),
  subtitle_text = "PC variance fraction x association effect size squared.",
  population = args$population
)
save_metric_barplot(
  plot_table = plot_selection_metrics,
  metric_column = "varpart_q3_percent",
  title_text = "variancePartition Q3 by covariate (75th percentile across genes)",
  y_label = "Q3 gene-level variance explained (%)",
  output_file = file.path(screen_plot_dir, "variance_partition_q3.png"),
  subtitle_text = "75th percentile across retained genes.",
  population = args$population
)
save_metric_barplot(
  plot_table = plot_selection_metrics,
  metric_column = "varpart_max_percent",
  title_text = "variancePartition max by covariate (largest single-gene value)",
  y_label = "Maximum gene-level variance explained (%)",
  output_file = file.path(screen_plot_dir, "variance_partition_max.png"),
  subtitle_text = "Largest value observed for any retained gene.",
  population = args$population
)
save_variance_partition_plot(
  varpart_object = variance_partition_results$varpart,
  output_file = file.path(screen_plot_dir, "variance_partition_multivariate.png"),
  population = args$population
)

save_selected_numeric_covariate_biology_plot(
  plot_data = numeric_biology_plot_data,
  biological_variable = "group_assignment",
  output_file = file.path(biology_balance_plot_dir, "selected_numeric_covariates_by_group_assignment.png"),
  population = args$population
)
save_selected_numeric_covariate_biology_plot(
  plot_data = numeric_biology_plot_data,
  biological_variable = "dpi",
  output_file = file.path(biology_balance_plot_dir, "selected_numeric_covariates_by_dpi.png"),
  population = args$population
)
save_selected_numeric_covariate_biology_plot(
  plot_data = numeric_biology_plot_data,
  biological_variable = "inoculum",
  output_file = file.path(biology_balance_plot_dir, "selected_numeric_covariates_by_inoculum.png"),
  population = args$population
)
save_selected_categorical_covariate_group_heatmaps(
  count_data = categorical_biology_plot_data,
  output_dir = biology_balance_plot_dir,
  population = args$population
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
    title_text = prefix_plot_title(args$population, sprintf("Blind PCA: %s", var_name)),
    output_file = file.path(screen_plot_dir, sprintf("blind_pca_by_%s.png", var_name))
  )
}

design_summary <- data.frame(
  design_id = "selected_covariates",
  design_formula = recommended_formula,
  stringsAsFactors = FALSE
)
write_tsv(design_summary, file.path(population_root, "designs_run.tsv"))

message("Finished covariate screen for population: ", args$population)
