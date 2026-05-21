#!/usr/bin/env Rscript

# Run covariate screening once across the full label-corrected dataset.
#
# This is a troubleshooting branch, not a replacement for the population-wise
# screen. The protected biological term is `population_group`, which preserves
# the population-specific group structure used for combined-model contrasts.

suppressPackageStartupMessages({
  required_packages <- c("dplyr", "DESeq2", "ggplot2", "pheatmap", "variancePartition", "here")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }
})

suppressMessages(here::i_am("scripts/troubleshooting/04h_whole_dataset_covariate_screening.R"))
source(here::here("scripts", "path_helpers.R"))
source(here::here("scripts", "plot_style.R"))
source(here::here("scripts", "02a_covariate_screening_utilities.R"))
source(here::here("scripts", "troubleshooting", "04a_deseq2_troubleshooting_utilities.R"))

parse_key_value_args_local <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %% 2 != 0) {
    stop("Arguments must be provided as --key value pairs.", call. = FALSE)
  }

  parsed <- defaults
  if (length(args) == 0) return(parsed)

  for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[[i]])
    value <- args[[i + 1]]
    if (!key %in% names(parsed)) {
      stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    }
    parsed[[key]] <- value
  }

  parsed
}

SETTINGS <- list(
  input_rds = trouble_input_rds_default,
  output_root = file.path(trouble_root_default, "07_whole_dataset_screened_model"),
  min_count = "10",
  min_samples = "auto",
  pc_variance_threshold = "0.85",
  pc_padj_cutoff = "0.05",
  weighted_pc_cutoff = "0.025",
  varpart_q3_cutoff = "0.03",
  varpart_max_cutoff = "0.75"
)

args <- parse_key_value_args_local(SETTINGS)
input_rds <- resolve_project_path(args$input_rds, must_work = TRUE)
output_root <- resolve_project_path(args$output_root)
min_count <- as.integer(args$min_count)
min_samples_arg <- args$min_samples

screen_root <- file.path(output_root, "covariate_screening")
screen_plot_dir <- file.path(screen_root, "plots")
biology_balance_plot_dir <- file.path(screen_plot_dir, "selected_covariates_by_biological_variables")
biology_balance_data_dir <- file.path(screen_root, "selected_covariates_by_biological_variables_data")
dir.create(screen_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(biology_balance_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(biology_balance_data_dir, recursive = TRUE, showWarnings = FALSE)

settings_numeric <- list(
  pc_variance_threshold = as.numeric(args$pc_variance_threshold),
  pc_padj_cutoff = as.numeric(args$pc_padj_cutoff),
  weighted_pc_cutoff = as.numeric(args$weighted_pc_cutoff),
  varpart_q3_cutoff = as.numeric(args$varpart_q3_cutoff),
  varpart_max_cutoff = as.numeric(args$varpart_max_cutoff)
)

analysis_input <- readRDS(input_rds)
required_objects <- c("txi_gene", "sample_metadata")
missing_objects <- setdiff(required_objects, names(analysis_input))
if (length(missing_objects) > 0) {
  stop("Input RDS missing object(s): ", paste(missing_objects, collapse = ", "), call. = FALSE)
}

txi_gene <- analysis_input$txi_gene
coldata <- prepare_coldata_for_modeling(analysis_input$sample_metadata)
coldata <- ensure_standardized_covariates(coldata)
coldata <- add_biological_fields(coldata)
rownames(coldata) <- coldata$sample

coldata$population_group <- factor(
  paste(as.character(coldata$population), as.character(coldata$group_assignment), sep = "__"),
  levels = as.vector(outer(exp383_population_levels, exp383_group_assignment_levels, paste, sep = "__"))
)
coldata$population_group <- droplevels(coldata$population_group)

min_samples <- if (identical(min_samples_arg, "auto")) {
  as.integer(min(table(coldata$population_group)))
} else {
  as.integer(min_samples_arg)
}

keep_genes <- compute_shared_filter(
  count_matrix = txi_gene$counts,
  grouping = coldata$population_group,
  min_count = min_count,
  min_samples = min_samples
)

txi_filtered <- subset_txi(
  txi = txi_gene,
  keep_rows = keep_genes,
  keep_cols = rep(TRUE, ncol(txi_gene$counts))
)

filter_summary <- data.frame(
  analysis_scope = "whole_dataset",
  n_samples = nrow(coldata),
  n_genes_before_filter = nrow(txi_gene$counts),
  n_genes_after_filter = nrow(txi_filtered$counts),
  min_count = min_count,
  min_samples = min_samples,
  filter_group_var = "population_group",
  stringsAsFactors = FALSE
)
write_tsv(filter_summary, file.path(output_root, "filter_summary.tsv"))

message("Running whole-dataset covariate screening")

registry <- build_covariate_registry(coldata)

# Replace the population-wise protected group term with the combined contrast
# term used by the downstream whole-dataset model.
registry <- registry[registry$covariate != "group_assignment", , drop = FALSE]
population_group_row <- data.frame(
  covariate = "population_group",
  analysis_column = "population_group",
  role = "protected_biological",
  is_candidate = FALSE,
  data_type = "categorical",
  n_non_missing = sum(!is.na(coldata$population_group)),
  n_unique = length(unique(stats::na.omit(coldata$population_group))),
  n_missing = sum(is.na(coldata$population_group)),
  stringsAsFactors = FALSE
)
registry <- rbind(population_group_row, registry)

pairwise_associations <- build_pairwise_covariate_associations(coldata, registry)
numeric_collinearity_summary <- build_numeric_collinearity_summary(pairwise_associations, registry)
formal_registry <- finalize_covariate_registry(registry, numeric_collinearity_summary, n_samples = nrow(coldata))

dds_blind <- DESeq2::DESeqDataSetFromTximport(
  txi = txi_filtered,
  colData = coldata,
  design = ~ 1
)
dds_blind <- DESeq2::estimateSizeFactors(dds_blind)
vsd_blind <- DESeq2::vst(dds_blind, blind = TRUE)

blind_pca_table <- compute_pca_table(vsd_blind, coldata)
percent_var <- attr(blind_pca_table, "percent_var")
selected_pcs <- select_pcs_by_cumulative_variance(
  percent_var,
  cumulative_threshold = settings_numeric$pc_variance_threshold
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
  pc_padj_cutoff = settings_numeric$pc_padj_cutoff
)

variance_partition_results <- run_multivariate_variance_partition(
  vsd = vsd_blind,
  registry = formal_registry,
  coldata = coldata
)
variance_partition_summary <- summarize_variance_partition(
  varpart_object = variance_partition_results$varpart,
  registry = formal_registry
)

selection_metrics <- apply_formal_selection_rule(
  registry = formal_registry,
  pc_summary = pc_covariate_summary,
  varpart_summary = variance_partition_summary,
  settings = settings_numeric
)
selected_covariates <- selection_metrics[selection_metrics$selected, , drop = FALSE]

selected_terms <- selected_covariates$analysis_column
recommended_formula <- paste("~", paste(unique(c(selected_terms, "population_group")), collapse = " + "))

biology_association_table <- build_selected_covariate_biology_associations(
  coldata = coldata,
  selection_metrics = selection_metrics
)
numeric_biology_plot_data <- build_selected_numeric_covariate_biology_data(
  coldata = coldata,
  selection_metrics = selection_metrics
)
categorical_biology_plot_data <- build_selected_categorical_covariate_biology_data(
  coldata = coldata,
  selection_metrics = selection_metrics
)

write_tsv(selection_metrics, file.path(screen_root, "covariate_selection_metrics.tsv"))
write_tsv(selected_covariates, file.path(screen_root, "selected_covariates.tsv"))
write_tsv(pairwise_associations, file.path(screen_root, "pairwise_covariate_associations.tsv"))
write_tsv(numeric_collinearity_summary, file.path(screen_root, "numeric_collinearity_summary.tsv"))
write_tsv(pc_variance_table, file.path(screen_root, "blind_pca_variance.tsv"))
write_tsv(pc_covariate_associations, file.path(screen_root, "pc_covariate_associations.tsv"))
write_tsv(pc_covariate_summary, file.path(screen_root, "pc_covariate_summary.tsv"))
write_tsv(variance_partition_summary, file.path(screen_root, "variance_partition_summary.tsv"))
write_tsv(variance_partition_results$input_summary, file.path(screen_root, "variance_partition_input_summary.tsv"))
write_tsv(biology_association_table, file.path(biology_balance_data_dir, "selected_covariate_biological_associations.tsv"))
write_tsv(numeric_biology_plot_data, file.path(biology_balance_data_dir, "selected_numeric_covariates_by_biological_variables.tsv"))
write_tsv(categorical_biology_plot_data, file.path(biology_balance_data_dir, "selected_categorical_covariates_by_biological_variables.tsv"))
write_tsv(
  data.frame(
    pc_variance_threshold = settings_numeric$pc_variance_threshold,
    pc_padj_cutoff = settings_numeric$pc_padj_cutoff,
    weighted_pc_cutoff = settings_numeric$weighted_pc_cutoff,
    varpart_q3_cutoff = settings_numeric$varpart_q3_cutoff,
    varpart_max_cutoff = settings_numeric$varpart_max_cutoff,
    stringsAsFactors = FALSE
  ),
  file.path(screen_root, "covariate_selection_thresholds.tsv")
)
writeLines(recommended_formula, con = file.path(screen_root, "recommended_design_formula.txt"))
write_tsv(
  data.frame(
    design_id = "whole_dataset_screened",
    design_formula = recommended_formula,
    stringsAsFactors = FALSE
  ),
  file.path(output_root, "designs_run.tsv")
)

save_pca_scree_plot(
  pc_variance_table = pc_variance_table,
  n_selected_pcs = length(selected_pcs),
  variance_threshold = settings_numeric$pc_variance_threshold,
  output_file = file.path(screen_plot_dir, "blind_pca_scree.png"),
  population = "All populations"
)
save_numeric_correlation_heatmap(
  coldata = coldata,
  registry = formal_registry,
  output_file = file.path(screen_plot_dir, "retained_numeric_covariate_correlations.png"),
  population = "All populations"
)

plot_selection_metrics <- selection_metrics[selection_metrics$is_candidate, , drop = FALSE]
save_metric_barplot(
  plot_table = plot_selection_metrics[!is.na(plot_selection_metrics$max_weighted_pc_score), , drop = FALSE],
  metric_column = "max_weighted_pc_score",
  title_text = "Maximum weighted PC association score",
  y_label = "Max weighted PC score",
  output_file = file.path(screen_plot_dir, "weighted_pc_scores.png"),
  subtitle_text = "PC variance fraction x association effect size squared.",
  population = "All populations"
)
save_metric_barplot(
  plot_table = plot_selection_metrics,
  metric_column = "varpart_q3_percent",
  title_text = "variancePartition Q3 by covariate (75th percentile across genes)",
  y_label = "Q3 gene-level variance explained (%)",
  output_file = file.path(screen_plot_dir, "variance_partition_q3.png"),
  subtitle_text = "75th percentile across retained genes.",
  population = "All populations"
)
save_metric_barplot(
  plot_table = plot_selection_metrics,
  metric_column = "varpart_max_percent",
  title_text = "variancePartition max by covariate (largest single-gene value)",
  y_label = "Maximum gene-level variance explained (%)",
  output_file = file.path(screen_plot_dir, "variance_partition_max.png"),
  subtitle_text = "Largest value observed for any retained gene.",
  population = "All populations"
)
save_variance_partition_plot(
  varpart_object = variance_partition_results$varpart,
  output_file = file.path(screen_plot_dir, "variance_partition_multivariate.png"),
  population = "All populations"
)

save_selected_numeric_covariate_biology_plot(
  plot_data = numeric_biology_plot_data,
  biological_variable = "group_assignment",
  output_file = file.path(biology_balance_plot_dir, "selected_numeric_covariates_by_group_assignment.png"),
  population = "All populations"
)
save_selected_numeric_covariate_biology_plot(
  plot_data = numeric_biology_plot_data,
  biological_variable = "dpi",
  output_file = file.path(biology_balance_plot_dir, "selected_numeric_covariates_by_dpi.png"),
  population = "All populations"
)
save_selected_numeric_covariate_biology_plot(
  plot_data = numeric_biology_plot_data,
  biological_variable = "inoculum",
  output_file = file.path(biology_balance_plot_dir, "selected_numeric_covariates_by_inoculum.png"),
  population = "All populations"
)
save_selected_categorical_covariate_group_heatmaps(
  count_data = categorical_biology_plot_data,
  output_dir = biology_balance_plot_dir,
  population = "All populations"
)

blind_pca_plot_table <- blind_pca_table
attr(blind_pca_plot_table, "percent_var") <- percent_var[1:2]
shared_plot_covariates <- unique(c(
  "population",
  "group_assignment",
  "population_group",
  selection_metrics$covariate[selection_metrics$selected]
))
for (var_name in shared_plot_covariates) {
  save_pca_plot(
    pca_table = blind_pca_plot_table,
    color_var = var_name,
    title_text = prefix_plot_title("All populations", sprintf("Blind PCA: %s", var_name)),
    output_file = file.path(screen_plot_dir, sprintf("blind_pca_by_%s.png", var_name))
  )
}

message("Finished whole-dataset covariate screening.")
