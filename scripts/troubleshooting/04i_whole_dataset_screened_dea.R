#!/usr/bin/env Rscript

# Fit and export a whole-dataset DESeq2 model using the whole-dataset screened
# technical covariates plus `population_group`.
#
# This keeps dispersion estimation shared across all samples while preserving
# population-specific contrasts through the population_group factor.

suppressPackageStartupMessages({
  required_packages <- c("dplyr", "DESeq2", "SummarizedExperiment", "ggplot2", "ggrepel", "pheatmap", "here")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }
})

suppressMessages(here::i_am("scripts/troubleshooting/04i_whole_dataset_screened_dea.R"))
source(here::here("scripts", "path_helpers.R"))
source(here::here("scripts", "plot_style.R"))
source(here::here("scripts", "03a_deseq2_dea_utilities.R"))
source(here::here("scripts", "troubleshooting", "04a_deseq2_troubleshooting_utilities.R"))

SETTINGS <- list(
  input_rds = trouble_input_rds_default,
  output_root = file.path(trouble_root_default, "07_whole_dataset_screened_model"),
  group_only_root = file.path(trouble_root_default, "04_combined_population_model"),
  contrast_tsv = trouble_contrast_tsv_default,
  min_count = as.character(trouble_min_count_default),
  min_samples = "auto",
  padj_cutoff = as.character(trouble_padj_cutoff_default),
  lfc_threshold = as.character(trouble_lfc_threshold_default),
  label_top_n = "20"
)

args <- parse_key_value_args(SETTINGS)
input_rds <- resolve_project_path(args$input_rds, must_work = TRUE)
output_root <- resolve_project_path(args$output_root)
group_only_root <- resolve_project_path(args$group_only_root, must_work = TRUE)
contrast_tsv <- resolve_project_path(args$contrast_tsv, must_work = TRUE)

screen_root <- file.path(output_root, "covariate_screening")
formula_file <- file.path(screen_root, "recommended_design_formula.txt")
stop_if_missing(formula_file, "Whole-dataset screened design formula")

design_formula <- trimws(readLines(formula_file, warn = FALSE))
design_formula <- design_formula[nzchar(design_formula)]
if (length(design_formula) != 1) {
  stop("Expected exactly one whole-dataset screened formula.", call. = FALSE)
}

min_count <- as.integer(args$min_count)
min_samples <- if (identical(args$min_samples, "auto")) "auto" else as.integer(args$min_samples)
padj_cutoff <- as.numeric(args$padj_cutoff)
lfc_threshold <- as.numeric(args$lfc_threshold)
label_top_n <- as.integer(args$label_top_n)

analysis_input <- load_dea_input(input_rds)
contrast_table <- read_contrast_table(contrast_tsv)
ensure_dir(output_root)

coldata <- add_biological_fields(prepare_coldata_for_deseq(analysis_input$sample_metadata))
coldata$population_group <- factor(
  paste(as.character(coldata$population), as.character(coldata$group_assignment), sep = "__"),
  levels = as.vector(outer(exp383_population_levels, exp383_group_assignment_levels, paste, sep = "__"))
)
coldata$population_group <- droplevels(coldata$population_group)
coldata <- ensure_formula_z_columns(coldata, design_formula)
rownames(coldata) <- coldata$sample

required_terms <- all.vars(stats::as.formula(design_formula))
missing_terms <- setdiff(required_terms, colnames(coldata))
if (length(missing_terms) > 0) {
  stop("Whole-dataset formula term(s) missing: ", paste(missing_terms, collapse = ", "), call. = FALSE)
}

min_samples_used <- if (identical(min_samples, "auto")) {
  as.integer(min(table(coldata$population_group)))
} else {
  as.integer(min_samples)
}

keep_genes <- compute_group_count_filter(
  count_matrix = analysis_input$txi_gene$counts,
  grouping = coldata$population_group,
  min_count = min_count,
  min_samples = min_samples_used
)
txi_filtered <- subset_txi(
  txi = analysis_input$txi_gene,
  keep_rows = keep_genes,
  keep_cols = rep(TRUE, ncol(analysis_input$txi_gene$counts))
)

rank_check <- check_design_matrix_full_rank(design_formula, coldata)
if (!rank_check$is_full_rank) {
  stop(
    sprintf(
      "Whole-dataset screened design matrix is not full rank: rank %s for %s columns.",
      rank_check$rank,
      rank_check$n_columns
    ),
    call. = FALSE
  )
}

model_root <- ensure_dir(file.path(output_root, "model_fit"))
plot_dir <- ensure_dir(file.path(model_root, "plots"))

message("Fitting whole-dataset screened DESeq2 model")
dds <- DESeq2::DESeqDataSetFromTximport(
  txi = txi_filtered,
  colData = coldata,
  design = stats::as.formula(design_formula)
)
dds <- DESeq2::DESeq(dds, quiet = TRUE)
vsd <- DESeq2::vst(dds, blind = FALSE)

saveRDS(dds, file.path(model_root, "dds.rds"))
saveRDS(vsd, file.path(model_root, "vst.rds"))
writeLines(design_formula, file.path(model_root, "model_formula.txt"))
writeLines(DESeq2::resultsNames(dds), file.path(model_root, "results_names.txt"))
write_tsv(
  data.frame(
    analysis_scope = "whole_dataset_screened",
    n_samples = ncol(dds),
    n_genes_before_filter = nrow(analysis_input$txi_gene$counts),
    n_genes_after_filter = nrow(dds),
    min_count = min_count,
    min_samples = min_samples_used,
    filter_group_var = "population_group",
    stringsAsFactors = FALSE
  ),
  file.path(model_root, "filter_summary.tsv")
)

design_matrix <- as.data.frame(rank_check$design_matrix, check.names = FALSE)
design_matrix$sample <- rownames(rank_check$design_matrix)
design_matrix <- design_matrix[, c("sample", setdiff(colnames(design_matrix), "sample")), drop = FALSE]
write_tsv(design_matrix, file.path(model_root, "design_matrix.tsv"))

sample_qc <- as.data.frame(SummarizedExperiment::colData(dds))
sample_qc$sample <- rownames(sample_qc)
sample_qc$size_factor <- get_sample_scaling_factors(dds)
sample_qc <- sample_qc[, c("sample", setdiff(colnames(sample_qc), "sample")), drop = FALSE]
write_tsv(sample_qc, file.path(model_root, "sample_model_qc.tsv"))

pca_table <- compute_pca_table(vsd, SummarizedExperiment::colData(dds))
write_tsv(pca_table, file.path(model_root, "vst_pca_scores.tsv"))
save_model_pca_plot(pca_table, "population", "Whole dataset screened model: VST PCA by population", file.path(plot_dir, "vst_pca_by_population.png"))
save_model_pca_plot(pca_table, "group_assignment", "Whole dataset screened model: VST PCA by group", file.path(plot_dir, "vst_pca_by_group_assignment.png"))
save_size_factor_plot(sample_qc, file.path(plot_dir, "size_factors.png"), "Whole dataset")
save_dispersion_plot(dds, file.path(plot_dir, "dispersion_estimates.png"), "Whole dataset")

model_summary <- data.frame(
  analysis_scope = "whole_dataset_screened",
  design_id = "whole_dataset_screened",
  n_samples = ncol(dds),
  n_genes = nrow(dds),
  design_formula = design_formula,
  design_rank = rank_check$rank,
  design_columns = rank_check$n_columns,
  min_count = min_count,
  min_samples = min_samples_used,
  stringsAsFactors = FALSE
)
write_tsv(model_summary, file.path(output_root, "model_fit_summary.tsv"))

baseline_table_path <- function(population, contrast_id) {
  file.path(
    group_only_root,
    "contrast_results",
    "tables",
    "all_genes",
    paste0(population, "_", contrast_id, ".tsv")
  )
}

read_result_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  result <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if ("gene_label" %in% colnames(result) && !"gene_symbol" %in% colnames(result)) {
    result$gene_symbol <- result$gene_label
  }
  result |>
    dplyr::arrange(is.na(.data$padj_zero), .data$padj_zero, .data$pvalue_zero)
}

summarise_result_hits <- function(result_table, prefix) {
  out <- data.frame(
    n_padj_lt_0_05 = sum(!is.na(result_table$padj_zero) & result_table$padj_zero < padj_cutoff),
    n_lfc_threshold_padj_lt_0_05 = sum(!is.na(result_table$padj_lfc_threshold) & result_table$padj_lfc_threshold < padj_cutoff),
    median_lfcSE = stats::median(result_table$lfcSE, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  names(out) <- paste0(prefix, "_", names(out))
  out
}

contrast_summaries <- list()
comparison_rows <- list()

for (population in exp383_population_levels) {
  population_root <- ensure_dir(file.path(output_root, population))
  result_root <- ensure_dir(file.path(population_root, "contrast_results"))
  all_gene_dir <- ensure_dir(file.path(result_root, "tables", "all_genes"))
  sig_zero_dir <- ensure_dir(file.path(result_root, "tables", "significant_genes"))
  sig_threshold_dir <- ensure_dir(file.path(result_root, "tables", "significant_lfc_threshold"))
  volcano_adjusted_dir <- ensure_dir(file.path(result_root, "plots", "volcano_adjusted_p"))
  volcano_raw_dir <- ensure_dir(file.path(result_root, "plots", "volcano_raw_p"))
  ma_dir <- ensure_dir(file.path(result_root, "plots", "ma"))

  for (i in seq_len(nrow(contrast_table))) {
    contrast_row <- contrast_table[i, , drop = FALSE]
    numerator <- paste(population, contrast_row$numerator_group, sep = "__")
    denominator <- paste(population, contrast_row$denominator_group, sep = "__")
    if (!all(c(numerator, denominator) %in% levels(dds$population_group))) next

    contrast_id <- contrast_row$contrast_id
    message("Exporting whole-dataset screened contrast ", population, " ", contrast_id)

    res_zero <- DESeq2::results(
      dds,
      contrast = c("population_group", numerator, denominator),
      alpha = padj_cutoff
    )
    res_threshold <- DESeq2::results(
      dds,
      contrast = c("population_group", numerator, denominator),
      alpha = padj_cutoff,
      lfcThreshold = lfc_threshold,
      altHypothesis = "greaterAbs"
    )

    result_table <- build_result_table(
      res_zero = res_zero,
      res_threshold = res_threshold,
      gene_annotation = analysis_input$gene_annotation,
      metadata_columns = list(
        population = population,
        design_id = "whole_dataset_screened",
        design_label = "whole-dataset screened covariates + population_group",
        contrast_id = contrast_id,
        contrast_family = contrast_row$contrast_family,
        contrast_label = contrast_row$contrast_label,
        numerator_group = numerator,
        denominator_group = denominator,
        model_formula = design_formula,
        lfc_threshold_test = lfc_threshold,
        lfc_threshold_alt_hypothesis = "greaterAbs"
      )
    )
    result_table$direction <- classify_de_direction(result_table$log2FoldChange, result_table$padj, padj_cutoff)
    result_table <- result_table |>
      dplyr::arrange(is.na(.data$padj_zero), .data$padj_zero, .data$pvalue_zero)

    sig_zero <- result_table[!is.na(result_table$padj_zero) & result_table$padj_zero < padj_cutoff, , drop = FALSE]
    sig_threshold <- result_table[result_table$direction != "not_significant" & !is.na(result_table$padj), , drop = FALSE]

    write_tsv(result_table, file.path(all_gene_dir, paste0(contrast_id, ".tsv")))
    write_tsv(sig_zero, file.path(sig_zero_dir, paste0(contrast_id, ".tsv")))
    write_tsv(sig_threshold, file.path(sig_threshold_dir, paste0(contrast_id, ".tsv")))

    plot_title <- paste0(population, ": ", contrast_row$contrast_label, " (whole-dataset screened)")
    save_volcano_plot(result_table, plot_title, file.path(volcano_adjusted_dir, paste0(contrast_id, ".png")), padj_cutoff, lfc_threshold, label_top_n, p_mode = "adjusted")
    save_volcano_plot(result_table, plot_title, file.path(volcano_raw_dir, paste0(contrast_id, ".png")), padj_cutoff, lfc_threshold, label_top_n, p_mode = "raw")
    save_ma_plot(result_table, plot_title, file.path(ma_dir, paste0(contrast_id, ".png")), padj_cutoff, lfc_threshold)

    contrast_summaries[[paste(population, contrast_id, sep = "__")]] <- summarise_result_table(
      result_table = result_table,
      population = population,
      contrast_row = contrast_row,
      design_id = "whole_dataset_screened",
      design_label = "whole-dataset screened covariates + population_group",
      model_formula = design_formula,
      padj_cutoff = padj_cutoff,
      lfc_threshold = lfc_threshold
    )

    baseline <- read_result_if_exists(baseline_table_path(population, contrast_id))
    if (!is.null(baseline)) {
      screened_zero_genes <- result_table$gene_id[!is.na(result_table$padj_zero) & result_table$padj_zero < padj_cutoff]
      baseline_zero_genes <- baseline$gene_id[!is.na(baseline$padj_zero) & baseline$padj_zero < padj_cutoff]
      screened_threshold_genes <- result_table$gene_id[!is.na(result_table$padj_lfc_threshold) & result_table$padj_lfc_threshold < padj_cutoff]
      baseline_threshold_genes <- baseline$gene_id[!is.na(baseline$padj_lfc_threshold) & baseline$padj_lfc_threshold < padj_cutoff]

      comparison_rows[[paste(population, contrast_id, sep = "__")]] <- cbind(
        data.frame(
          population = population,
          contrast = contrast_id,
          stringsAsFactors = FALSE
        ),
        summarise_result_hits(result_table, "screened"),
        summarise_result_hits(baseline, "group_only"),
        data.frame(
          zero_effect_overlap = length(intersect(screened_zero_genes, baseline_zero_genes)),
          zero_effect_screened_only = length(setdiff(screened_zero_genes, baseline_zero_genes)),
          zero_effect_group_only_only = length(setdiff(baseline_zero_genes, screened_zero_genes)),
          lfc_threshold_overlap = length(intersect(screened_threshold_genes, baseline_threshold_genes)),
          lfc_threshold_screened_only = length(setdiff(screened_threshold_genes, baseline_threshold_genes)),
          lfc_threshold_group_only_only = length(setdiff(baseline_threshold_genes, screened_threshold_genes)),
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

contrast_summary <- do.call(rbind, contrast_summaries)
rownames(contrast_summary) <- NULL
write_tsv(contrast_summary, file.path(output_root, "whole_dataset_screened_contrast_summary.tsv"))

if (length(comparison_rows) > 0) {
  comparison <- do.call(rbind, comparison_rows)
  rownames(comparison) <- NULL
  comparison$delta_zero_effect_n_padj_lt_0_05 <- comparison$screened_n_padj_lt_0_05 - comparison$group_only_n_padj_lt_0_05
  comparison$delta_lfc_threshold_n_padj_lt_0_05 <- comparison$screened_n_lfc_threshold_padj_lt_0_05 - comparison$group_only_n_lfc_threshold_padj_lt_0_05
  comparison$delta_median_lfcSE <- comparison$screened_median_lfcSE - comparison$group_only_median_lfcSE
  write_tsv(comparison, file.path(output_root, "whole_dataset_screened_vs_group_only_comparison.tsv"))
}

combine_significant_tables <- function(subdir, output_name) {
  files <- Sys.glob(file.path(output_root, "*", "contrast_results", "tables", subdir, "*.tsv"))
  tables <- lapply(files, function(path) {
    result <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(result) == 0) return(NULL)
    result$source_file <- path
    result
  })
  tables <- Filter(Negate(is.null), tables)
  combined <- if (length(tables) == 0) {
    data.frame()
  } else {
    dplyr::bind_rows(tables) |>
      dplyr::arrange(is.na(.data$padj_zero), .data$padj_zero, .data$pvalue_zero, .data$population, .data$contrast_id)
  }

  write_tsv(combined, file.path(output_root, output_name))
}

combine_significant_tables("significant_lfc_threshold", "whole_dataset_screened_all_significant_lfc_threshold.tsv")
combine_significant_tables("significant_genes", "whole_dataset_screened_all_significant_zero_effect.tsv")

message("Finished whole-dataset screened DEA.")
