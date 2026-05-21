#!/usr/bin/env Rscript

# Export population-wise DESeq2 DEA results for configured group contrasts.
#
# Inputs:
# - fitted DESeq2 objects from 03b_deseq2_dea_fit_models.R
# - explicit contrast definitions from config/dea_contrasts.tsv
#
# Outputs:
# - all-gene DESeq2 result tables, including zero-effect and threshold-test p-values
# - significant-gene tables called by the threshold test
# - adjusted-p and raw-p volcano plots with ggrepel labels
# - MA plots
# - compact contrast summary tables

suppressPackageStartupMessages({
  required_packages <- c("DESeq2", "SummarizedExperiment", "dplyr", "ggplot2", "ggrepel", "here")
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

suppressMessages(here::i_am("scripts/03c_deseq2_dea_export_results.R"))
source(here::here("scripts", "path_helpers.R"))
source(here::here("scripts", "plot_style.R"))
source(here::here("scripts", "03a_deseq2_dea_utilities.R"))

# ----------------
# Fixed settings
# ----------------
# By default, output is written beside the fitted models produced by 03b.
SETTINGS <- list(
  population = "all",
  input_rds = "results/dea/01e_sample_label_correction/exp383_salmon_gene_input_label_corrected.rds",
  fit_root = "results/dea/03_deseq2_dea",
  output_root = "results/dea/03_deseq2_dea",
  contrast_tsv = "config/dea_contrasts.tsv",
  padj_cutoff = "0.05",
  lfc_threshold = "0.5",
  label_top_n = "20"
)

args <- parse_key_value_args(SETTINGS)

input_rds <- resolve_project_path(args$input_rds, must_work = TRUE)
fit_root <- resolve_project_path(args$fit_root, must_work = TRUE)
output_root <- resolve_project_path(args$output_root)
contrast_tsv <- resolve_project_path(args$contrast_tsv, must_work = TRUE)

padj_cutoff <- as.numeric(args$padj_cutoff)
lfc_threshold <- as.numeric(args$lfc_threshold)
label_top_n <- as.integer(args$label_top_n)

if (is.na(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff >= 1) {
  stop("--padj_cutoff must be between 0 and 1.", call. = FALSE)
}
if (is.na(lfc_threshold) || lfc_threshold < 0) {
  stop("--lfc_threshold must be a non-negative number.", call. = FALSE)
}
if (is.na(label_top_n) || label_top_n < 0) {
  stop("--label_top_n must be a non-negative integer.", call. = FALSE)
}

# ----------------
# Load inputs
# ----------------
analysis_input <- load_dea_input(input_rds)
contrast_table <- read_contrast_table(contrast_tsv)

available_populations <- list.dirs(fit_root, full.names = FALSE, recursive = FALSE)
available_populations <- intersect(order_known_levels(available_populations, exp383_population_levels), available_populations)
populations_to_run <- parse_population_selection(args$population, available_populations)

ensure_dir(output_root)

# ---------------------------
# Export configured contrasts
# ---------------------------
combined_summaries <- list()

for (population in populations_to_run) {
  message("Exporting DESeq2 contrasts for ", population)

  model_root <- file.path(fit_root, population, "model_fit")
  dds_file <- file.path(model_root, "dds.rds")
  stop_if_missing(dds_file, sprintf("%s fitted DESeq2 object", population))

  dds <- readRDS(dds_file)
  group_levels <- levels(droplevels(SummarizedExperiment::colData(dds)$group_assignment))

  runnable <- contrast_table$numerator_group %in% group_levels &
    contrast_table$denominator_group %in% group_levels
  population_contrasts <- contrast_table[runnable, , drop = FALSE]
  skipped_contrasts <- contrast_table[!runnable, , drop = FALSE]

  if (nrow(population_contrasts) == 0) {
    warning(sprintf("No configured contrasts are runnable for %s.", population), call. = FALSE)
    next
  }

  result_root <- ensure_dir(file.path(output_root, population, "contrast_results"))
  all_gene_dir <- ensure_dir(file.path(result_root, "tables", "all_genes"))
  significant_gene_dir <- ensure_dir(file.path(result_root, "tables", "significant_lfc_threshold"))
  volcano_adjusted_dir <- ensure_dir(file.path(result_root, "plots", "volcano_adjusted_p"))
  volcano_raw_dir <- ensure_dir(file.path(result_root, "plots", "volcano_raw_p"))
  ma_dir <- ensure_dir(file.path(result_root, "plots", "ma"))

  if (nrow(skipped_contrasts) > 0) {
    write_tsv(skipped_contrasts, file.path(result_root, "skipped_contrasts.tsv"))
  }

  population_summaries <- list()

  for (i in seq_len(nrow(population_contrasts))) {
    contrast_row <- population_contrasts[i, , drop = FALSE]
    contrast_id <- contrast_row$contrast_id

    # DESeq2 contrast direction is numerator minus denominator.
    # Keep the conventional zero-effect test for auditability.
    res_zero <- DESeq2::results(
      dds,
      contrast = c(
        "group_assignment",
        contrast_row$numerator_group,
        contrast_row$denominator_group
      ),
      alpha = padj_cutoff
    )

    # Primary calling uses a non-zero threshold test. This asks whether the
    # absolute effect is significantly larger than the configured log2FC value,
    # rather than testing against zero and filtering after the fact.
    res_threshold <- DESeq2::results(
      dds,
      contrast = c(
        "group_assignment",
        contrast_row$numerator_group,
        contrast_row$denominator_group
      ),
      alpha = padj_cutoff,
      lfcThreshold = lfc_threshold,
      altHypothesis = "greaterAbs"
    )

    zero_table <- as.data.frame(res_zero)
    threshold_table <- as.data.frame(res_threshold)
    zero_table$gene_id <- rownames(zero_table)
    threshold_table$gene_id <- rownames(threshold_table)

    result_table <- data.frame(
      gene_id = zero_table$gene_id,
      baseMean = zero_table$baseMean,
      log2FoldChange = zero_table$log2FoldChange,
      lfcSE = zero_table$lfcSE,
      stat_zero = zero_table$stat,
      pvalue_zero = zero_table$pvalue,
      padj_zero = zero_table$padj,
      stat_lfc_threshold = threshold_table$stat,
      pvalue_lfc_threshold = threshold_table$pvalue,
      padj_lfc_threshold = threshold_table$padj,
      stringsAsFactors = FALSE
    )

    # Downstream plotting functions use pvalue/padj as the active test result.
    result_table$pvalue <- result_table$pvalue_lfc_threshold
    result_table$padj <- result_table$padj_lfc_threshold
    rownames(result_table) <- result_table$gene_id

    result_table <- annotate_deseq_results(
      result_table = result_table,
      gene_annotation = analysis_input$gene_annotation
    )
    result_table$population <- population
    result_table$contrast_id <- contrast_id
    result_table$contrast_family <- contrast_row$contrast_family
    result_table$contrast_label <- contrast_row$contrast_label
    result_table$numerator_group <- contrast_row$numerator_group
    result_table$denominator_group <- contrast_row$denominator_group
    result_table$lfc_threshold_test <- lfc_threshold
    result_table$lfc_threshold_alt_hypothesis <- "greaterAbs"
    result_table$direction <- classify_de_direction(
      log2_fold_change = result_table$log2FoldChange,
      padj = result_table$padj,
      padj_cutoff = padj_cutoff
    )

    result_table <- result_table |>
      dplyr::arrange(is.na(.data$padj_zero), .data$padj_zero, .data$pvalue_zero)

    significant_table <- result_table[
      result_table$direction != "not_significant" &
        !is.na(result_table$padj),
      ,
      drop = FALSE
    ]

    write_tsv(result_table, file.path(all_gene_dir, paste0(contrast_id, ".tsv")))
    write_tsv(significant_table, file.path(significant_gene_dir, paste0(contrast_id, ".tsv")))

    plot_title <- paste0(population, ": ", contrast_row$contrast_label)

    save_volcano_plot(
      result_table = result_table,
      contrast_label = plot_title,
      output_file = file.path(volcano_adjusted_dir, paste0(contrast_id, ".png")),
      padj_cutoff = padj_cutoff,
      lfc_threshold = lfc_threshold,
      label_top_n = label_top_n,
      p_mode = "adjusted"
    )

    save_volcano_plot(
      result_table = result_table,
      contrast_label = plot_title,
      output_file = file.path(volcano_raw_dir, paste0(contrast_id, ".png")),
      padj_cutoff = padj_cutoff,
      lfc_threshold = lfc_threshold,
      label_top_n = label_top_n,
      p_mode = "raw"
    )

    save_ma_plot(
      result_table = result_table,
      contrast_label = plot_title,
      output_file = file.path(ma_dir, paste0(contrast_id, ".png")),
      padj_cutoff = padj_cutoff,
      lfc_threshold = lfc_threshold
    )

    population_summaries[[contrast_id]] <- data.frame(
      population = population,
      contrast_id = contrast_id,
      contrast_family = contrast_row$contrast_family,
      contrast_label = contrast_row$contrast_label,
      numerator_group = contrast_row$numerator_group,
      denominator_group = contrast_row$denominator_group,
      n_genes_total = nrow(result_table),
      n_genes_tested_zero = sum(!is.na(result_table$padj_zero)),
      n_sig_zero_padj = sum(!is.na(result_table$padj_zero) & result_table$padj_zero <= padj_cutoff),
      n_genes_tested_lfc_threshold = sum(!is.na(result_table$padj_lfc_threshold)),
      n_sig_lfc_threshold_padj = sum(!is.na(result_table$padj_lfc_threshold) & result_table$padj_lfc_threshold <= padj_cutoff),
      n_sig_padj = sum(!is.na(result_table$padj_lfc_threshold) & result_table$padj_lfc_threshold <= padj_cutoff),
      n_sig_up = sum(result_table$direction == "up"),
      n_sig_down = sum(result_table$direction == "down"),
      padj_cutoff = padj_cutoff,
      lfc_threshold = lfc_threshold,
      lfc_threshold_alt_hypothesis = "greaterAbs",
      stringsAsFactors = FALSE
    )
  }

  population_summary <- do.call(rbind, population_summaries)
  rownames(population_summary) <- NULL
  write_tsv(population_summary, file.path(result_root, "contrast_summary.tsv"))
  combined_summaries[[population]] <- population_summary
}

if (length(combined_summaries) > 0) {
  combined_summary <- do.call(rbind, combined_summaries)
  rownames(combined_summary) <- NULL
  write_tsv(combined_summary, file.path(output_root, "contrast_summary.tsv"))
}

message("Finished DESeq2 contrast export.")
