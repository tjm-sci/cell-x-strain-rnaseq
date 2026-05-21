#!/usr/bin/env Rscript

# Build objective QC exclusion sets and rerun selected DESeq2 sensitivity fits.

suppressPackageStartupMessages({
  required_packages <- c("dplyr", "DESeq2", "SummarizedExperiment", "ggplot2", "ggrepel", "pheatmap", "here")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }
})

suppressMessages(here::i_am("scripts/troubleshooting/04d_sample_qc_sensitivity.R"))
source(here::here("scripts", "path_helpers.R"))
source(here::here("scripts", "plot_style.R"))
source(here::here("scripts", "03a_deseq2_dea_utilities.R"))
source(here::here("scripts", "troubleshooting", "04a_deseq2_troubleshooting_utilities.R"))

SETTINGS <- list(
  population = "all",
  input_rds = trouble_input_rds_default,
  screening_root = trouble_screening_root_default,
  primary_fit_root = "results/dea/03_deseq2_dea",
  confounding_root = file.path(trouble_root_default, "02_covariate_confounding"),
  output_root = file.path(trouble_root_default, "03_sample_qc_sensitivity"),
  contrast_tsv = trouble_contrast_tsv_default,
  min_count = as.character(trouble_min_count_default),
  min_samples = "auto",
  padj_cutoff = as.character(trouble_padj_cutoff_default),
  lfc_threshold = as.character(trouble_lfc_threshold_default),
  label_top_n = "20"
)

args <- parse_key_value_args(SETTINGS)

input_rds <- resolve_project_path(args$input_rds, must_work = TRUE)
screening_root <- resolve_project_path(args$screening_root, must_work = TRUE)
primary_fit_root <- resolve_project_path(args$primary_fit_root, must_work = TRUE)
confounding_root <- resolve_project_path(args$confounding_root, must_work = TRUE)
output_root <- resolve_project_path(args$output_root)
contrast_tsv <- resolve_project_path(args$contrast_tsv, must_work = TRUE)

min_count <- as.integer(args$min_count)
min_samples <- if (identical(args$min_samples, "auto")) "auto" else as.integer(args$min_samples)
padj_cutoff <- as.numeric(args$padj_cutoff)
lfc_threshold <- as.numeric(args$lfc_threshold)
label_top_n <- as.integer(args$label_top_n)

analysis_input <- load_dea_input(input_rds)
metadata <- add_biological_fields(prepare_coldata_for_deseq(analysis_input$sample_metadata))
available_populations <- order_known_levels(metadata$population, exp383_population_levels)
populations_to_run <- parse_population_selection(args$population, available_populations)
contrast_table <- read_contrast_table(contrast_tsv)

ensure_dir(output_root)
exclusion_root <- ensure_dir(file.path(output_root, "exclusion_sets"))

missing_qc_columns <- character()
threshold_rows <- list()
audit_rows <- list()
exclusion_sets <- list()

safe_read_tsv <- function(path) {
  stop_if_missing(path, path)
  utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

compute_robust_pc_outliers <- function(pca_table, pc_cols) {
  pc_matrix <- as.matrix(pca_table[, pc_cols, drop = FALSE])
  center <- apply(pc_matrix, 2, stats::median, na.rm = TRUE)
  scale <- apply(pc_matrix, 2, stats::mad, na.rm = TRUE)
  scale[is.na(scale) | scale == 0] <- 1
  robust_distance <- sqrt(rowSums(sweep(sweep(pc_matrix, 2, center, "-"), 2, scale, "/")^2))
  out <- iqr_outlier_flags(robust_distance, low = FALSE, high = TRUE)
  list(distance = robust_distance, flags = out$flags, lower = out$lower, upper = out$upper)
}

safe_safer_formula <- function(population, coldata) {
  risk_file <- file.path(confounding_root, "covariate_risk_classification.tsv")
  risk_table <- safe_read_tsv(risk_file)
  risk_pop <- risk_table[risk_table$population == population, , drop = FALSE]
  keep <- risk_pop$risk_class %in% c("low_risk", "moderate_risk")
  terms <- risk_pop$analysis_column[keep]
  term_available <- function(term) {
    term %in% colnames(coldata) ||
      (grepl("_z$", term) && sub("_z$", "", term) %in% colnames(coldata))
  }
  terms <- terms[vapply(terms, term_available, logical(1))]
  formula_from_terms(c(terms, "group_assignment"))
}

for (population in populations_to_run) {
  message("Building sample QC audit for ", population)
  sample_qc <- safe_read_tsv(file.path(primary_fit_root, population, "model_fit", "sample_model_qc.tsv"))
  pca_table <- safe_read_tsv(file.path(primary_fit_root, population, "model_fit", "vst_pca_scores.tsv"))
  vsd <- readRDS(file.path(primary_fit_root, population, "model_fit", "vst.rds"))

  coldata <- droplevels(metadata[metadata$population == population, , drop = FALSE])
  rownames(coldata) <- coldata$sample
  counts_pop <- analysis_input$txi_gene$counts[, coldata$sample, drop = FALSE]
  n_expressed <- colSums(counts_pop >= min_count)

  audit <- data.frame(
    sample_id = coldata$sample,
    animal_id = coldata$mouse_n,
    population = population,
    group_id = as.character(coldata$group_assignment),
    inoculum = as.character(coldata$inoculum),
    dpi = as.character(coldata$dpi),
    infected_status = as.character(coldata$infected_status),
    stringsAsFactors = FALSE
  )

  qc_columns <- c(
    "salmon_percent_mapped",
    "salmon_num_mapped_millions",
    "rrna_alignment_percent",
    "filtered_total_sequences_millions"
  )
  for (col in qc_columns) {
    if (col %in% colnames(coldata)) {
      audit[[col]] <- as.numeric(coldata[[col]])
    } else {
      audit[[col]] <- NA_real_
      missing_qc_columns <- c(missing_qc_columns, paste(population, col, sep = ":"))
    }
  }
  audit$n_expressed_genes <- as.integer(n_expressed[audit$sample_id])
  audit$size_factor <- sample_qc$size_factor[match(audit$sample_id, sample_qc$sample)]
  audit$pca1 <- pca_table$PC1[match(audit$sample_id, pca_table$sample)]
  audit$pca2 <- pca_table$PC2[match(audit$sample_id, pca_table$sample)]

  pc_cols <- grep("^PC[0-9]+$", colnames(pca_table), value = TRUE)
  pc_cols <- pc_cols[seq_len(min(length(pc_cols), 5L))]
  pca_out <- compute_robust_pc_outliers(pca_table[match(audit$sample_id, pca_table$sample), , drop = FALSE], pc_cols)

  assay_matrix <- SummarizedExperiment::assay(vsd)
  sample_dist <- as.matrix(stats::dist(t(assay_matrix)))
  mean_distance <- rowMeans(sample_dist, na.rm = TRUE)
  mean_distance <- mean_distance[audit$sample_id]
  distance_out <- iqr_outlier_flags(mean_distance, low = FALSE, high = TRUE)

  salmon_low <- iqr_outlier_flags(audit$salmon_percent_mapped, low = TRUE, high = FALSE)
  rrna_high <- iqr_outlier_flags(audit$rrna_alignment_percent, low = FALSE, high = TRUE)
  size_out <- iqr_outlier_flags(audit$size_factor, low = TRUE, high = TRUE)

  audit$pca_robust_distance <- pca_out$distance
  audit$mean_sample_distance <- mean_distance
  audit$pca_outlier_flag <- pca_out$flags
  audit$sample_distance_outlier_flag <- distance_out$flags
  audit$qc_outlier_flag <- salmon_low$flags | rrna_high$flags | size_out$flags

  reasons <- vector("list", nrow(audit))
  for (i in seq_len(nrow(audit))) {
    r <- character()
    if (salmon_low$flags[[i]]) r <- c(r, "low_salmon_percent_mapped")
    if (rrna_high$flags[[i]]) r <- c(r, "high_rrna_alignment_percent")
    if (size_out$flags[[i]]) r <- c(r, "extreme_size_factor")
    if (audit$pca_outlier_flag[[i]]) r <- c(r, "pca_robust_distance_outlier")
    if (audit$sample_distance_outlier_flag[[i]]) r <- c(r, "sample_distance_outlier")
    reasons[[i]] <- paste(r, collapse = ";")
  }
  audit$outlier_reason <- unlist(reasons)

  audit_rows[[population]] <- audit

  threshold_rows[[paste(population, "salmon", sep = "__")]] <- data.frame(population = population, metric = "salmon_percent_mapped", direction = "low", lower_threshold = salmon_low$lower, upper_threshold = salmon_low$upper, multiplier = 3, stringsAsFactors = FALSE)
  threshold_rows[[paste(population, "rrna", sep = "__")]] <- data.frame(population = population, metric = "rrna_alignment_percent", direction = "high", lower_threshold = rrna_high$lower, upper_threshold = rrna_high$upper, multiplier = 3, stringsAsFactors = FALSE)
  threshold_rows[[paste(population, "size", sep = "__")]] <- data.frame(population = population, metric = "size_factor", direction = "both", lower_threshold = size_out$lower, upper_threshold = size_out$upper, multiplier = 3, stringsAsFactors = FALSE)
  threshold_rows[[paste(population, "pca", sep = "__")]] <- data.frame(population = population, metric = "pca_robust_distance_PC1_to_PC5", direction = "high", lower_threshold = pca_out$lower, upper_threshold = pca_out$upper, multiplier = 3, stringsAsFactors = FALSE)
  threshold_rows[[paste(population, "sample_distance", sep = "__")]] <- data.frame(population = population, metric = "mean_sample_distance", direction = "high", lower_threshold = distance_out$lower, upper_threshold = distance_out$upper, multiplier = 3, stringsAsFactors = FALSE)

  exclusion_sets[[population]] <- list(
    set_00_no_exclusion = character(),
    set_01_extreme_qc_outliers_only = audit$sample_id[audit$qc_outlier_flag],
    set_02_pca_distance_outliers_only = audit$sample_id[audit$pca_outlier_flag | audit$sample_distance_outlier_flag],
    set_03_combined_extreme_qc_and_pca_outliers = audit$sample_id[audit$qc_outlier_flag | audit$pca_outlier_flag | audit$sample_distance_outlier_flag]
  )

  population_exclusion_table <- do.call(
    rbind,
    lapply(names(exclusion_sets[[population]]), function(set_id) {
      samples <- exclusion_sets[[population]][[set_id]]
      data.frame(
        population = population,
        exclusion_set = set_id,
        n_samples_removed = length(samples),
        samples_removed = paste(samples, collapse = ";"),
        stringsAsFactors = FALSE
      )
    })
  )
  write_tsv(population_exclusion_table, file.path(exclusion_root, paste0(population, "_exclusion_sets.tsv")))
}

audit_all <- do.call(rbind, audit_rows)
rownames(audit_all) <- NULL
write_tsv(audit_all, file.path(output_root, "sample_qc_audit_all_populations.tsv"))

thresholds <- do.call(rbind, threshold_rows)
rownames(thresholds) <- NULL
write_tsv(thresholds, file.path(exclusion_root, "exclusion_thresholds.tsv"))

if (length(missing_qc_columns) == 0) {
  writeLines("No requested QC columns were missing.", file.path(output_root, "missing_qc_columns.txt"))
} else {
  writeLines(unique(missing_qc_columns), file.path(output_root, "missing_qc_columns.txt"))
}

all_sensitivity_summaries <- list()
skipped_fits <- list()

for (population in populations_to_run) {
  message("Running QC sensitivity fits for ", population)
  coldata_pop <- droplevels(metadata[metadata$population == population, , drop = FALSE])
  original_formula <- read_population_design_formula(screening_root, population)
  safer_formula <- safe_safer_formula(population, coldata_pop)

  fit_grid <- data.frame(
    exclusion_set = c(
      "set_00_no_exclusion",
      "set_00_no_exclusion",
      "set_03_combined_extreme_qc_and_pca_outliers",
      "set_03_combined_extreme_qc_and_pca_outliers"
    ),
    design_id = c(
      "minimal_group",
      "original_screened_design",
      "minimal_group",
      "safer_selected_design"
    ),
    design_label = c(
      "minimal group",
      "original screened design",
      "minimal group after combined QC/PCA exclusion",
      "safer selected design after combined QC/PCA exclusion"
    ),
    model_formula = c("~ group_assignment", original_formula, "~ group_assignment", safer_formula),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(fit_grid))) {
    row <- fit_grid[i, , drop = FALSE]
    exclude_samples <- exclusion_sets[[population]][[row$exclusion_set]]
    n_removed <- length(exclude_samples)
    remaining <- coldata_pop[!coldata_pop$sample %in% exclude_samples, , drop = FALSE]
    remaining_group_n <- table(droplevels(remaining$group_assignment))

    notes <- "fit completed"
    if (any(remaining_group_n < 2)) {
      notes <- "skipped: exclusion set leaves a group with fewer than two samples"
      skipped_fits[[paste(population, row$exclusion_set, row$design_id, sep = "__")]] <- data.frame(
        population = population,
        exclusion_set = row$exclusion_set,
        design_id = row$design_id,
        notes = notes,
        stringsAsFactors = FALSE
      )
      next
    }

    fit_root <- ensure_dir(file.path(output_root, population, row$exclusion_set, row$design_id))
    result <- tryCatch(
      {
        prepared <- prepare_population_troubleshooting(
          analysis_input = analysis_input,
          population = population,
          design_formula = row$model_formula,
          min_count = min_count,
          min_samples = min_samples,
          exclude_samples = exclude_samples
        )
        fit_troubleshooting_model(
          prepared = prepared,
          design_formula = row$model_formula,
          population = population,
          output_dir = fit_root,
          design_id = row$design_id,
          design_label = row$design_label
        )
        dds <- readRDS(file.path(fit_root, "model_fit", "dds.rds"))
        contrast_summary <- export_troubleshooting_contrasts(
          dds = dds,
          analysis_input = analysis_input,
          contrast_table = contrast_table,
          output_dir = fit_root,
          population = population,
          design_id = row$design_id,
          design_label = row$design_label,
          model_formula = row$model_formula,
          padj_cutoff = padj_cutoff,
          lfc_threshold = lfc_threshold,
          label_top_n = label_top_n
        )
        contrast_summary
      },
      error = function(e) {
        notes <<- paste("skipped:", conditionMessage(e))
        NULL
      }
    )

    if (is.null(result)) {
      skipped_fits[[paste(population, row$exclusion_set, row$design_id, sep = "__")]] <- data.frame(
        population = population,
        exclusion_set = row$exclusion_set,
        design_id = row$design_id,
        notes = notes,
        stringsAsFactors = FALSE
      )
      next
    }

    result$exclusion_set <- row$exclusion_set
    result$n_samples_removed <- n_removed
    result$samples_removed <- paste(exclude_samples, collapse = ";")
    result$notes <- notes
    result <- result[, c(
      "population", "exclusion_set", "n_samples_removed", "samples_removed",
      "design_id", "design_label", "model_formula", "contrast", "contrast_label",
      "n_genes_tested", "n_nominal_p_lt_0_05", "n_padj_lt_0_1", "n_padj_lt_0_05",
      "n_lfc_threshold_padj_lt_0_05", "median_lfcSE", "notes",
      setdiff(colnames(result), c(
        "population", "exclusion_set", "n_samples_removed", "samples_removed",
        "design_id", "design_label", "model_formula", "contrast", "contrast_label",
        "n_genes_tested", "n_nominal_p_lt_0_05", "n_padj_lt_0_1", "n_padj_lt_0_05",
        "n_lfc_threshold_padj_lt_0_05", "median_lfcSE", "notes"
      ))
    )]
    all_sensitivity_summaries[[paste(population, row$exclusion_set, row$design_id, sep = "__")]] <- result
  }
}

if (length(all_sensitivity_summaries) > 0) {
  sensitivity_summary <- do.call(rbind, all_sensitivity_summaries)
  rownames(sensitivity_summary) <- NULL
  write_tsv(sensitivity_summary, file.path(output_root, "sample_qc_sensitivity_contrast_summary.tsv"))
}

if (length(skipped_fits) > 0) {
  skipped <- do.call(rbind, skipped_fits)
  rownames(skipped) <- NULL
  write_tsv(skipped, file.path(output_root, "skipped_sensitivity_fits.tsv"))
}

message("Finished sample QC sensitivity analysis.")
