#!/usr/bin/env Rscript

# Audit expected positive-control genes across PU1 design-ablation outputs.

suppressPackageStartupMessages({
  required_packages <- c("dplyr", "DESeq2", "SummarizedExperiment", "ggplot2", "here")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }
})

suppressMessages(here::i_am("scripts/troubleshooting/04f_positive_control_gene_audit.R"))
source(here::here("scripts", "path_helpers.R"))
source(here::here("scripts", "plot_style.R"))
source(here::here("scripts", "03a_deseq2_dea_utilities.R"))
source(here::here("scripts", "troubleshooting", "04a_deseq2_troubleshooting_utilities.R"))

SETTINGS <- list(
  input_rds = trouble_input_rds_default,
  design_ablation_root = file.path(trouble_root_default, "01_deseq2_design_ablation"),
  output_root = file.path(trouble_root_default, "05_positive_control_gene_audit"),
  gene_panel_tsv = "config/troubleshooting_positive_control_genes.tsv",
  population = "PU1",
  padj_cutoff = as.character(trouble_padj_cutoff_default),
  nominal_cutoff = "0.05",
  lfc_small_threshold = "0.5"
)

args <- parse_key_value_args(SETTINGS)
input_rds <- resolve_project_path(args$input_rds, must_work = TRUE)
design_ablation_root <- resolve_project_path(args$design_ablation_root, must_work = TRUE)
output_root <- resolve_project_path(args$output_root)
gene_panel_tsv <- resolve_project_path(args$gene_panel_tsv, must_work = TRUE)
population <- args$population
padj_cutoff <- as.numeric(args$padj_cutoff)
nominal_cutoff <- as.numeric(args$nominal_cutoff)
lfc_small_threshold <- as.numeric(args$lfc_small_threshold)

analysis_input <- load_dea_input(input_rds)
gene_panel <- utils::read.delim(gene_panel_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
gene_panel <- gene_panel[gene_panel$expected_population == population | is.na(gene_panel$expected_population), , drop = FALSE]
gene_annotation <- as.data.frame(analysis_input$gene_annotation, stringsAsFactors = FALSE)

ensure_dir(output_root)
plot_root <- ensure_dir(file.path(output_root, "plots", population))

panel_map <- merge(
  gene_panel,
  gene_annotation,
  by.x = "gene_symbol",
  by.y = "gene_name",
  all.x = TRUE,
  sort = FALSE
)

design_table <- utils::read.delim(
  file.path(design_ablation_root, population, "design_ladder.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

classify_positive_control <- function(row, expected_direction, detected, passed_prefilter) {
  if (is.na(row$gene_id) || row$gene_id == "") return("missing_gene_mapping")
  if (!detected) return("not_detected")
  if (!passed_prefilter) return("filtered_out")
  if (!is.na(row$log2FoldChange)) {
    if (identical(expected_direction, "up") && row$log2FoldChange < -0.25) return("opposite_direction")
    if (identical(expected_direction, "down") && row$log2FoldChange > 0.25) return("opposite_direction")
  }
  if (!is.na(row$padj) && row$padj < padj_cutoff) return("fdr_significant")
  if (!is.na(row$pvalue) && row$pvalue < nominal_cutoff) return("nominal_only")
  if (!is.na(row$baseMean) && row$baseMean < 10) return("detected_low_baseMean")
  if (!is.na(row$log2FoldChange) && abs(row$log2FoldChange) < lfc_small_threshold) return("small_lfc")
  if (!is.na(row$lfcSE) && (row$lfcSE >= 1 || (!is.na(row$log2FoldChange) && row$lfcSE >= abs(row$log2FoldChange)))) {
    return("large_lfc_high_se")
  }
  "large_lfc_high_se"
}

all_rows <- list()
count_matrix <- analysis_input$txi_gene$counts

for (i in seq_len(nrow(design_table))) {
  design_id <- design_table$design_id[[i]]
  design_label <- design_table$design_label[[i]]
  all_gene_dir <- file.path(design_ablation_root, population, design_id, "contrast_results", "tables", "all_genes")
  if (!dir.exists(all_gene_dir)) next
  contrast_files <- list.files(all_gene_dir, pattern = "\\.tsv$", full.names = TRUE)

  for (contrast_file in contrast_files) {
    contrast_id <- sub("\\.tsv$", "", basename(contrast_file))
    result_table <- utils::read.delim(contrast_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

    for (j in seq_len(nrow(panel_map))) {
      gene_symbol <- panel_map$gene_symbol[[j]]
      gene_id <- panel_map$gene_id[[j]]
      detected <- !is.na(gene_id) && gene_id %in% rownames(count_matrix) && sum(count_matrix[gene_id, ], na.rm = TRUE) > 0
      passed <- !is.na(gene_id) && gene_id %in% result_table$gene_id
      result_row <- if (passed) {
        result_table[result_table$gene_id == gene_id, , drop = FALSE][1, , drop = FALSE]
      } else {
        data.frame(
          gene_id = gene_id,
          baseMean = NA_real_,
          log2FoldChange = NA_real_,
          lfcSE = NA_real_,
          stat_zero = NA_real_,
          pvalue_zero = NA_real_,
          padj_zero = NA_real_,
          stringsAsFactors = FALSE
        )
      }

      audit_row <- data.frame(
        population = population,
        design_id = design_id,
        design_label = design_label,
        contrast = contrast_id,
        gene_symbol = gene_symbol,
        gene_id = gene_id,
        baseMean = result_row$baseMean,
        log2FoldChange = result_row$log2FoldChange,
        lfcSE = result_row$lfcSE,
        stat = result_row$stat_zero,
        pvalue = result_row$pvalue_zero,
        padj = result_row$padj_zero,
        padj_is_na = is.na(result_row$padj_zero),
        pvalue_is_na = is.na(result_row$pvalue_zero),
        detected_in_count_matrix = detected,
        passed_prefilter = passed,
        expected_direction = panel_map$expected_direction[[j]],
        classification = NA_character_,
        stringsAsFactors = FALSE
      )
      audit_row$classification <- classify_positive_control(
        audit_row,
        expected_direction = panel_map$expected_direction[[j]],
        detected = detected,
        passed_prefilter = passed
      )
      all_rows[[paste(design_id, contrast_id, gene_symbol, gene_id, sep = "__")]] <- audit_row
    }
  }
}

audit <- do.call(rbind, all_rows)
rownames(audit) <- NULL
write_tsv(audit, file.path(output_root, paste0(population, "_positive_control_gene_audit.tsv")))

summary_rows <- list()
for (key in unique(paste(audit$population, audit$design_id, audit$contrast, sep = "__"))) {
  idx <- paste(audit$population, audit$design_id, audit$contrast, sep = "__") == key
  x <- audit[idx, , drop = FALSE]
  summary_rows[[key]] <- data.frame(
    population = x$population[[1]],
    design_id = x$design_id[[1]],
    contrast = x$contrast[[1]],
    n_positive_control_genes_in_matrix = sum(x$detected_in_count_matrix, na.rm = TRUE),
    n_not_detected = sum(x$classification == "not_detected", na.rm = TRUE),
    n_filtered_out = sum(x$classification == "filtered_out", na.rm = TRUE),
    n_nominal_only = sum(x$classification == "nominal_only", na.rm = TRUE),
    n_fdr_significant = sum(x$classification == "fdr_significant", na.rm = TRUE),
    n_large_lfc_high_se = sum(x$classification == "large_lfc_high_se", na.rm = TRUE),
    n_opposite_direction = sum(x$classification == "opposite_direction", na.rm = TRUE),
    median_positive_control_lfc = stats::median(x$log2FoldChange, na.rm = TRUE),
    median_positive_control_lfcSE = stats::median(x$lfcSE, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
summary_table <- do.call(rbind, summary_rows)
rownames(summary_table) <- NULL
write_tsv(summary_table, file.path(output_root, "positive_control_summary_by_design.tsv"))

# Plot normalised counts from the PU1 minimal design. This keeps the visual audit
# independent of the covariate-adjusted design choices being tested.
minimal_dds_file <- file.path(design_ablation_root, population, "design_01_minimal_group", "model_fit", "dds.rds")
if (file.exists(minimal_dds_file)) {
  dds <- readRDS(minimal_dds_file)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  coldata <- as.data.frame(SummarizedExperiment::colData(dds))
  coldata$sample <- rownames(coldata)

  for (gene_symbol in unique(panel_map$gene_symbol)) {
    gene_ids <- unique(panel_map$gene_id[panel_map$gene_symbol == gene_symbol])
    gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids %in% rownames(norm_counts)]
    if (length(gene_ids) == 0) next
    count_values <- colSums(norm_counts[gene_ids, , drop = FALSE])
    plot_df <- data.frame(
      sample = names(count_values),
      normalised_count = as.numeric(count_values),
      group_assignment = coldata$group_assignment[match(names(count_values), coldata$sample)],
      stringsAsFactors = FALSE
    )
    plot_df$group_assignment <- factor(as.character(plot_df$group_assignment), levels = exp383_group_assignment_levels)
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = group_assignment, y = log10(normalised_count + 1), fill = group_assignment)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.75, show.legend = FALSE) +
      ggplot2::geom_jitter(width = 0.16, height = 0, size = 1.8, alpha = 0.85, show.legend = FALSE) +
      exp383_theme(base_size = 11) +
      ggplot2::scale_fill_manual(values = exp383_group_assignment_palette, labels = format_group_assignment_label) +
      ggplot2::labs(
        title = paste0(population, ": ", gene_symbol),
        subtitle = "DESeq2 normalised counts from the minimal PU1 design.",
        x = "Group",
        y = "log10(normalised count + 1)"
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
    exp383_save_ggplot(file.path(plot_root, paste0(sanitize_filename(gene_symbol), "_normalised_counts_by_group.png")), plot = p, width = 8, height = 5.5)
  }
}

message("Finished positive-control gene audit.")
