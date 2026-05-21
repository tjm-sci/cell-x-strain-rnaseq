#!/usr/bin/env Rscript

# Fit a combined population_group DESeq2 sensitivity model.

suppressPackageStartupMessages({
  required_packages <- c("dplyr", "DESeq2", "SummarizedExperiment", "ggplot2", "here")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }
})

suppressMessages(here::i_am("scripts/troubleshooting/04e_combined_population_model.R"))
source(here::here("scripts", "path_helpers.R"))
source(here::here("scripts", "plot_style.R"))
source(here::here("scripts", "03a_deseq2_dea_utilities.R"))
source(here::here("scripts", "troubleshooting", "04a_deseq2_troubleshooting_utilities.R"))

SETTINGS <- list(
  input_rds = trouble_input_rds_default,
  design_ablation_root = file.path(trouble_root_default, "01_deseq2_design_ablation"),
  primary_dea_root = "results/dea/03_deseq2_dea",
  output_root = file.path(trouble_root_default, "04_combined_population_model"),
  contrast_tsv = trouble_contrast_tsv_default,
  min_count = as.character(trouble_min_count_default),
  min_samples = "auto",
  padj_cutoff = as.character(trouble_padj_cutoff_default),
  lfc_threshold = as.character(trouble_lfc_threshold_default)
)

args <- parse_key_value_args(SETTINGS)
input_rds <- resolve_project_path(args$input_rds, must_work = TRUE)
design_ablation_root <- resolve_project_path(args$design_ablation_root)
primary_dea_root <- resolve_project_path(args$primary_dea_root, must_work = TRUE)
output_root <- resolve_project_path(args$output_root)
contrast_tsv <- resolve_project_path(args$contrast_tsv, must_work = TRUE)

min_count <- as.integer(args$min_count)
min_samples <- if (identical(args$min_samples, "auto")) "auto" else as.integer(args$min_samples)
padj_cutoff <- as.numeric(args$padj_cutoff)
lfc_threshold <- as.numeric(args$lfc_threshold)

analysis_input <- load_dea_input(input_rds)
contrast_table <- read_contrast_table(contrast_tsv)
ensure_dir(output_root)

coldata <- add_biological_fields(prepare_coldata_for_deseq(analysis_input$sample_metadata))
rownames(coldata) <- coldata$sample
coldata$population_group <- factor(
  paste(as.character(coldata$population), as.character(coldata$group_assignment), sep = "__"),
  levels = as.vector(outer(exp383_population_levels, exp383_group_assignment_levels, paste, sep = "__"))
)
coldata$population_group <- droplevels(coldata$population_group)

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

model_root <- ensure_dir(file.path(output_root, "model_fit"))
message("Fitting combined population_group model")
dds <- DESeq2::DESeqDataSetFromTximport(
  txi = txi_filtered,
  colData = coldata,
  design = ~ population_group
)
dds <- DESeq2::DESeq(dds, quiet = TRUE)
vsd <- DESeq2::vst(dds, blind = FALSE)
saveRDS(dds, file.path(model_root, "dds.rds"))
saveRDS(vsd, file.path(model_root, "vst.rds"))
writeLines("~ population_group", file.path(model_root, "model_formula.txt"))
writeLines(DESeq2::resultsNames(dds), file.path(model_root, "results_names.txt"))
write_tsv(
  data.frame(
    population = "combined",
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

sample_qc <- as.data.frame(SummarizedExperiment::colData(dds))
sample_qc$sample <- rownames(sample_qc)
sample_qc$size_factor <- get_sample_scaling_factors(dds)
sample_qc <- sample_qc[, c("sample", setdiff(colnames(sample_qc), "sample")), drop = FALSE]
write_tsv(sample_qc, file.path(model_root, "sample_model_qc.tsv"))
pca_table <- compute_pca_table(vsd, SummarizedExperiment::colData(dds))
write_tsv(pca_table, file.path(model_root, "vst_pca_scores.tsv"))

result_root <- ensure_dir(file.path(output_root, "contrast_results"))
all_gene_dir <- ensure_dir(file.path(result_root, "tables", "all_genes"))

split_original_summary <- utils::read.delim(file.path(primary_dea_root, "contrast_summary.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
split_minimal_summary_file <- file.path(design_ablation_root, "design_ablation_contrast_summary.tsv")
split_minimal_summary <- if (file.exists(split_minimal_summary_file)) {
  x <- utils::read.delim(split_minimal_summary_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  x[x$design_id == "design_01_minimal_group", , drop = FALSE]
} else {
  data.frame()
}

primary_all_gene_path <- function(population, contrast_id) {
  file.path(primary_dea_root, population, "contrast_results", "tables", "all_genes", paste0(contrast_id, ".tsv"))
}

median_lfcse_from_table <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  x <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  stats::median(x$lfcSE, na.rm = TRUE)
}

combined_summaries <- list()
for (population in exp383_population_levels) {
  for (i in seq_len(nrow(contrast_table))) {
    contrast_row <- contrast_table[i, , drop = FALSE]
    numerator <- paste(population, contrast_row$numerator_group, sep = "__")
    denominator <- paste(population, contrast_row$denominator_group, sep = "__")
    if (!all(c(numerator, denominator) %in% levels(dds$population_group))) next

    contrast_id <- paste(population, contrast_row$contrast_id, sep = "_")
    message("Exporting combined contrast ", contrast_id)

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
        contrast_id = contrast_id,
        contrast_label = paste0(population, ": ", contrast_row$contrast_label),
        numerator_group = numerator,
        denominator_group = denominator,
        model_formula = "~ population_group",
        lfc_threshold_test = lfc_threshold,
        lfc_threshold_alt_hypothesis = "greaterAbs"
      )
    )
    result_table$direction <- classify_de_direction(result_table$log2FoldChange, result_table$padj, padj_cutoff)
    result_table <- result_table[order(is.na(result_table$padj_zero), result_table$padj_zero, result_table$pvalue_zero), , drop = FALSE]
    write_tsv(result_table, file.path(all_gene_dir, paste0(contrast_id, ".tsv")))

    split_original_row <- split_original_summary[
      split_original_summary$population == population & split_original_summary$contrast_id == contrast_row$contrast_id,
      ,
      drop = FALSE
    ]
    split_minimal_row <- split_minimal_summary[
      split_minimal_summary$population == population & split_minimal_summary$contrast == contrast_row$contrast_id,
      ,
      drop = FALSE
    ]

    combined_summaries[[contrast_id]] <- data.frame(
      population = population,
      contrast = contrast_row$contrast_id,
      combined_model_n_padj_lt_0_05 = sum(!is.na(result_table$padj_zero) & result_table$padj_zero < 0.05),
      split_original_n_padj_lt_0_05 = if (nrow(split_original_row) == 1) split_original_row$n_sig_zero_padj else NA_integer_,
      split_minimal_n_padj_lt_0_05 = if (nrow(split_minimal_row) == 1) split_minimal_row$n_padj_lt_0_05 else NA_integer_,
      combined_model_median_lfcSE = stats::median(result_table$lfcSE, na.rm = TRUE),
      split_original_median_lfcSE = median_lfcse_from_table(primary_all_gene_path(population, contrast_row$contrast_id)),
      split_minimal_median_lfcSE = if (nrow(split_minimal_row) == 1) split_minimal_row$median_lfcSE else NA_real_,
      notes = "Combined sensitivity model uses ~ population_group.",
      stringsAsFactors = FALSE
    )
  }
}

summary_table <- do.call(rbind, combined_summaries)
rownames(summary_table) <- NULL
write_tsv(summary_table, file.path(output_root, "combined_population_contrast_summary.tsv"))

message("Finished combined population model.")
