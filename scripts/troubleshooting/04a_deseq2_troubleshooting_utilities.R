# Shared helpers for EXP383 DESeq2 troubleshooting scripts.
#
# These functions intentionally build on the main DEA utilities rather than
# creating a second pipeline. All outputs are written under
# results/dea/troubleshooting/ so primary DEA outputs are preserved.

trouble_root_default <- "results/dea/troubleshooting"
trouble_input_rds_default <- "results/dea/01e_sample_label_correction/exp383_salmon_gene_input_label_corrected.rds"
trouble_screening_root_default <- "results/dea/02_covariate_screening"
trouble_contrast_tsv_default <- "config/dea_contrasts.tsv"
trouble_min_count_default <- 10L
trouble_padj_cutoff_default <- 0.05
trouble_lfc_threshold_default <- 0.5

trouble_qc_covariates <- c(
  "salmon_percent_mapped",
  "salmon_num_mapped_millions",
  "filtered_total_sequences_millions",
  "fastp_passed_filter_reads_millions",
  "raw_total_sequences_millions",
  "rrna_alignment_percent",
  "fastp_percent_adapter",
  "fastp_percent_duplication",
  "filtered_percent_gc",
  "raw_percent_gc",
  "filtered_percent_duplicates",
  "raw_percent_duplicates",
  "fastp_q30_rate_after_filtering",
  "fastp_percent_surviving"
)

trouble_direct_depth_mapping_covariates <- c(
  "salmon_percent_mapped",
  "salmon_num_mapped_millions",
  "filtered_total_sequences_millions",
  "raw_total_sequences_millions",
  "fastp_passed_filter_reads_millions",
  "fastp_percent_surviving"
)

trouble_biological_vars <- c("group_assignment", "inoculum", "dpi", "infected_status")

add_biological_fields <- function(coldata) {
  coldata <- as.data.frame(coldata, stringsAsFactors = FALSE)
  if ("group_assignment" %in% colnames(coldata)) {
    coldata$group_assignment <- factor(
      as.character(coldata$group_assignment),
      levels = order_known_levels(coldata$group_assignment, exp383_group_assignment_levels)
    )
    coldata$group_id <- coldata$group_assignment
  }
  if ("inoculum" %in% colnames(coldata)) {
    coldata$inoculum <- factor(as.character(coldata$inoculum), levels = exp383_inoculum_levels)
  }
  if ("dpi" %in% colnames(coldata)) {
    coldata$dpi <- factor(as.character(coldata$dpi), levels = exp383_dpi_levels)
  }
  if ("population" %in% colnames(coldata)) {
    coldata$population <- factor(as.character(coldata$population), levels = exp383_population_levels)
  }
  if ("inoculum" %in% colnames(coldata)) {
    coldata$infected_status <- ifelse(as.character(coldata$inoculum) == "CBH", "CBH", "infected")
    coldata$infected_status <- factor(coldata$infected_status, levels = c("CBH", "infected"))
  }
  coldata
}

formula_from_terms <- function(terms) {
  terms <- terms[nzchar(terms)]
  terms <- unique(terms)
  paste("~", paste(terms, collapse = " + "))
}

read_selected_covariates <- function(screening_root, population) {
  cov_file <- file.path(screening_root, population, "covariate_screening", "selected_covariates.tsv")
  stop_if_missing(cov_file, sprintf("%s selected covariates", population))
  utils::read.delim(cov_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

selected_analysis_columns <- function(selected_covariates, roles = NULL) {
  selected <- selected_covariates
  if ("selected" %in% colnames(selected)) {
    selected <- selected[selected$selected %in% c(TRUE, "TRUE", "true", "1", 1), , drop = FALSE]
  }
  if (!is.null(roles) && "covariate_role" %in% colnames(selected)) {
    selected <- selected[selected$covariate_role %in% roles, , drop = FALSE]
  }
  cols <- selected$analysis_column
  cols[is.na(cols) | !nzchar(cols)] <- selected$covariate[is.na(cols) | !nzchar(cols)]
  unique(cols[nzchar(cols)])
}

selected_non_depth_mapping_terms <- function(selected_covariates, coldata) {
  selected <- selected_covariates
  if ("selected" %in% colnames(selected)) {
    selected <- selected[selected$selected %in% c(TRUE, "TRUE", "true", "1", 1), , drop = FALSE]
  }

  terms <- selected$analysis_column
  terms[is.na(terms) | !nzchar(terms)] <- selected$covariate[is.na(terms) | !nzchar(terms)]

  source_terms <- sub("_z$", "", terms)
  keep <- !source_terms %in% trouble_direct_depth_mapping_covariates
  terms <- unique(terms[keep & nzchar(terms)])

  term_available <- function(term) {
    term %in% colnames(coldata) ||
      (grepl("_z$", term) && sub("_z$", "", term) %in% colnames(coldata))
  }
  terms[vapply(terms, term_available, logical(1))]
}

build_design_ladder <- function(screening_root, population, coldata) {
  selected <- read_selected_covariates(screening_root, population)
  original_formula <- read_population_design_formula(screening_root, population)

  date_cols <- intersect(c("date_nuc_prep_days_z", "date_nuc_prep_days"), colnames(coldata))
  date_selected <- selected$covariate %in% "date_nuc_prep_days" | selected$analysis_column %in% date_cols

  term_available <- function(term) {
    term %in% colnames(coldata) ||
      (grepl("_z$", term) && sub("_z$", "", term) %in% colnames(coldata))
  }

  non_qc_terms <- selected_analysis_columns(selected, roles = "technical_metadata")
  non_qc_terms <- non_qc_terms[vapply(non_qc_terms, term_available, logical(1))]
  qc_terms <- selected_analysis_columns(selected, roles = "technical_qc")
  qc_terms <- qc_terms[vapply(qc_terms, term_available, logical(1))]
  non_depth_mapping_terms <- selected_non_depth_mapping_terms(selected, coldata)

  designs <- list(
    data.frame(
      design_id = "design_01_minimal_group",
      design_label = "minimal group",
      model_formula = "~ group_assignment",
      stringsAsFactors = FALSE
    ),
    data.frame(
      design_id = "design_02_known_batch_group",
      design_label = "inoculation batch + group",
      model_formula = "~ inoculation_batch + group_assignment",
      stringsAsFactors = FALSE
    )
  )

  if (any(date_selected) && length(date_cols) > 0) {
    designs[[length(designs) + 1L]] <- data.frame(
      design_id = "design_03_date_group",
      design_label = "nuclei prep date + group",
      model_formula = formula_from_terms(c(date_cols[[1]], "group_assignment")),
      stringsAsFactors = FALSE
    )
  }

  if (length(non_qc_terms) > 0) {
    designs[[length(designs) + 1L]] <- data.frame(
      design_id = "design_04_selected_non_qc_covariates_group",
      design_label = "selected non-QC covariates + group",
      model_formula = formula_from_terms(c(non_qc_terms, "group_assignment")),
      stringsAsFactors = FALSE
    )
  }

  if (length(qc_terms) > 0) {
    designs[[length(designs) + 1L]] <- data.frame(
      design_id = "design_05_selected_qc_covariates_group",
      design_label = "selected QC covariates + group",
      model_formula = formula_from_terms(c(qc_terms, "group_assignment")),
      stringsAsFactors = FALSE
    )
  }

  designs[[length(designs) + 1L]] <- data.frame(
    design_id = "design_06_original_screened_design",
    design_label = "original screened design",
    model_formula = original_formula,
    stringsAsFactors = FALSE
  )

  non_depth_mapping_formula <- formula_from_terms(c(non_depth_mapping_terms, "group_assignment"))
  existing_formulas <- vapply(designs, function(x) x$model_formula[[1]], character(1))
  if (length(non_depth_mapping_terms) > 0 && !non_depth_mapping_formula %in% existing_formulas) {
    designs[[length(designs) + 1L]] <- data.frame(
      design_id = "design_07_selected_non_depth_mapping_covariates_group",
      design_label = "selected non-depth/non-Salmon-mapping covariates + group",
      model_formula = non_depth_mapping_formula,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, designs)
  out <- out[!duplicated(out[, c("design_id", "model_formula")]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

prepare_population_troubleshooting <- function(analysis_input, population, design_formula,
                                               min_count = trouble_min_count_default,
                                               min_samples = "auto",
                                               exclude_samples = character()) {
  coldata_all <- prepare_coldata_for_deseq(analysis_input$sample_metadata)
  coldata_all <- add_biological_fields(coldata_all)
  coldata_all <- ensure_formula_z_columns(coldata_all, design_formula)

  keep_samples <- as.character(coldata_all$population) == population
  if (length(exclude_samples) > 0) {
    keep_samples <- keep_samples & !as.character(coldata_all$sample) %in% exclude_samples
  }
  coldata_pop <- droplevels(coldata_all[keep_samples, , drop = FALSE])
  rownames(coldata_pop) <- coldata_pop$sample

  required_terms <- all.vars(stats::as.formula(design_formula))
  missing_terms <- setdiff(required_terms, colnames(coldata_pop))
  if (length(missing_terms) > 0) {
    stop(sprintf("%s design term(s) missing: %s", population, paste(missing_terms, collapse = ", ")), call. = FALSE)
  }

  if (!"group_assignment" %in% colnames(coldata_pop)) {
    stop("group_assignment is required for troubleshooting contrasts.", call. = FALSE)
  }

  txi_pop <- subset_txi(
    txi = analysis_input$txi_gene,
    keep_rows = rep(TRUE, nrow(analysis_input$txi_gene$counts)),
    keep_cols = keep_samples
  )

  min_samples_used <- if (identical(min_samples, "auto")) {
    as.integer(min(table(coldata_pop$group_assignment)))
  } else {
    as.integer(min_samples)
  }

  keep_genes <- compute_group_count_filter(
    count_matrix = txi_pop$counts,
    grouping = coldata_pop$group_assignment,
    min_count = min_count,
    min_samples = min_samples_used
  )

  txi_filtered <- subset_txi(
    txi = txi_pop,
    keep_rows = keep_genes,
    keep_cols = rep(TRUE, ncol(txi_pop$counts))
  )

  filter_summary <- data.frame(
    population = population,
    n_samples = nrow(coldata_pop),
    n_genes_before_filter = nrow(txi_pop$counts),
    n_genes_after_filter = nrow(txi_filtered$counts),
    min_count = min_count,
    min_samples = min_samples_used,
    filter_group_var = "group_assignment",
    n_excluded_samples = length(exclude_samples),
    excluded_samples = paste(exclude_samples, collapse = ";"),
    stringsAsFactors = FALSE
  )

  list(coldata = coldata_pop, txi = txi_filtered, keep_genes = keep_genes, filter_summary = filter_summary)
}

fit_troubleshooting_model <- function(prepared, design_formula, population, output_dir,
                                      design_id = NA_character_, design_label = NA_character_) {
  rank_check <- check_design_matrix_full_rank(design_formula, prepared$coldata)
  if (!rank_check$is_full_rank) {
    stop(
      sprintf(
        "%s %s design matrix is not full rank: rank %s for %s columns.",
        population,
        design_id,
        rank_check$rank,
        rank_check$n_columns
      ),
      call. = FALSE
    )
  }

  model_root <- ensure_dir(file.path(output_dir, "model_fit"))
  plot_dir <- ensure_dir(file.path(model_root, "plots"))

  dds <- DESeq2::DESeqDataSetFromTximport(
    txi = prepared$txi,
    colData = prepared$coldata,
    design = stats::as.formula(design_formula)
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  vsd <- DESeq2::vst(dds, blind = FALSE)

  saveRDS(dds, file.path(model_root, "dds.rds"))
  saveRDS(vsd, file.path(model_root, "vst.rds"))
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

  pca_table <- compute_pca_table(vsd, SummarizedExperiment::colData(dds))
  write_tsv(pca_table, file.path(model_root, "vst_pca_scores.tsv"))

  save_model_pca_plot(pca_table, "group_assignment", paste0(population, ": ", design_label, " model PCA"), file.path(plot_dir, "vst_pca_group_assignment.png"))
  save_size_factor_plot(sample_qc, file.path(plot_dir, "size_factors.png"), population)
  save_dispersion_plot(dds, file.path(plot_dir, "dispersion_estimates.png"), population)

  data.frame(
    population = population,
    design_id = design_id,
    design_label = design_label,
    n_samples = ncol(dds),
    n_genes = nrow(dds),
    design_formula = design_formula,
    design_rank = rank_check$rank,
    design_columns = rank_check$n_columns,
    min_count = prepared$filter_summary$min_count,
    min_samples = prepared$filter_summary$min_samples,
    stringsAsFactors = FALSE
  )
}

build_result_table <- function(res_zero, res_threshold, gene_annotation, metadata_columns = list()) {
  zero_table <- as.data.frame(res_zero, stringsAsFactors = FALSE, check.names = FALSE)
  threshold_table <- as.data.frame(res_threshold, stringsAsFactors = FALSE, check.names = FALSE)
  zero_table$gene_id <- rownames(zero_table)
  threshold_table$gene_id <- rownames(threshold_table)

  result_table <- dplyr::tibble(
    gene_id = zero_table$gene_id,
    baseMean = zero_table$baseMean,
    log2FoldChange = zero_table$log2FoldChange,
    lfcSE = zero_table$lfcSE,
    stat_zero = zero_table$stat,
    pvalue_zero = zero_table$pvalue,
    padj_zero = zero_table$padj,
    stat_lfc_threshold = threshold_table$stat,
    pvalue_lfc_threshold = threshold_table$pvalue,
    padj_lfc_threshold = threshold_table$padj
  )
  result_table$pvalue <- result_table$pvalue_lfc_threshold
  result_table$padj <- result_table$padj_lfc_threshold

  result_table <- annotate_deseq_results(result_table, gene_annotation)
  for (col_name in names(metadata_columns)) {
    result_table[[col_name]] <- metadata_columns[[col_name]]
  }

  result_table |>
    dplyr::arrange(is.na(.data$padj_zero), .data$padj_zero, .data$pvalue_zero) |>
    dplyr::select("gene_symbol", "gene_id", dplyr::everything())
}

summarise_result_table <- function(result_table, population, contrast_row, design_id, design_label,
                                   model_formula, padj_cutoff, lfc_threshold) {
  data.frame(
    population = population,
    design_id = design_id,
    design_label = design_label,
    model_formula = model_formula,
    contrast = contrast_row$contrast_id,
    contrast_label = contrast_row$contrast_label,
    numerator_group = contrast_row$numerator_group,
    denominator_group = contrast_row$denominator_group,
    n_genes_tested = sum(!is.na(result_table$pvalue_zero)),
    n_nominal_p_lt_0_05 = sum(!is.na(result_table$pvalue_zero) & result_table$pvalue_zero < 0.05),
    n_padj_lt_0_1 = sum(!is.na(result_table$padj_zero) & result_table$padj_zero < 0.1),
    n_padj_lt_0_05 = sum(!is.na(result_table$padj_zero) & result_table$padj_zero < 0.05),
    n_padj_lt_0_01 = sum(!is.na(result_table$padj_zero) & result_table$padj_zero < 0.01),
    n_lfc_threshold_padj_lt_0_05 = sum(!is.na(result_table$padj_lfc_threshold) & result_table$padj_lfc_threshold < padj_cutoff),
    median_abs_lfc = stats::median(abs(result_table$log2FoldChange), na.rm = TRUE),
    median_lfcSE = stats::median(result_table$lfcSE, na.rm = TRUE),
    n_genes_with_padj_NA = sum(is.na(result_table$padj_zero)),
    n_genes_with_pvalue_NA = sum(is.na(result_table$pvalue_zero)),
    padj_cutoff = padj_cutoff,
    lfc_threshold = lfc_threshold,
    stringsAsFactors = FALSE
  )
}

export_troubleshooting_contrasts <- function(dds, analysis_input, contrast_table, output_dir,
                                             population, design_id, design_label, model_formula,
                                             padj_cutoff = trouble_padj_cutoff_default,
                                             lfc_threshold = trouble_lfc_threshold_default,
                                             label_top_n = 20L) {
  result_root <- ensure_dir(file.path(output_dir, "contrast_results"))
  all_gene_dir <- ensure_dir(file.path(result_root, "tables", "all_genes"))
  sig_zero_dir <- ensure_dir(file.path(result_root, "tables", "significant_genes"))
  sig_threshold_dir <- ensure_dir(file.path(result_root, "tables", "significant_lfc_threshold"))
  volcano_adjusted_dir <- ensure_dir(file.path(result_root, "plots", "volcano_adjusted_p"))
  volcano_raw_dir <- ensure_dir(file.path(result_root, "plots", "volcano_raw_p"))
  ma_dir <- ensure_dir(file.path(result_root, "plots", "ma"))

  group_levels <- levels(droplevels(SummarizedExperiment::colData(dds)$group_assignment))
  runnable <- contrast_table$numerator_group %in% group_levels &
    contrast_table$denominator_group %in% group_levels
  population_contrasts <- contrast_table[runnable, , drop = FALSE]
  skipped_contrasts <- contrast_table[!runnable, , drop = FALSE]
  if (nrow(skipped_contrasts) > 0) {
    write_tsv(skipped_contrasts, file.path(result_root, "skipped_contrasts.tsv"))
  }

  summaries <- list()
  for (i in seq_len(nrow(population_contrasts))) {
    contrast_row <- population_contrasts[i, , drop = FALSE]
    res_zero <- DESeq2::results(
      dds,
      contrast = c("group_assignment", contrast_row$numerator_group, contrast_row$denominator_group),
      alpha = padj_cutoff
    )
    res_threshold <- DESeq2::results(
      dds,
      contrast = c("group_assignment", contrast_row$numerator_group, contrast_row$denominator_group),
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
        design_id = design_id,
        design_label = design_label,
        model_formula = model_formula,
        contrast_id = contrast_row$contrast_id,
        contrast_family = contrast_row$contrast_family,
        contrast_label = contrast_row$contrast_label,
        numerator_group = contrast_row$numerator_group,
        denominator_group = contrast_row$denominator_group,
        lfc_threshold_test = lfc_threshold,
        lfc_threshold_alt_hypothesis = "greaterAbs"
      )
    )
    result_table$direction <- classify_de_direction(result_table$log2FoldChange, result_table$padj, padj_cutoff)
    result_table <- result_table |>
      dplyr::arrange(is.na(.data$padj_zero), .data$padj_zero, .data$pvalue_zero)

    sig_zero <- result_table[!is.na(result_table$padj_zero) & result_table$padj_zero < padj_cutoff, , drop = FALSE]
    sig_threshold <- result_table[result_table$direction != "not_significant" & !is.na(result_table$padj), , drop = FALSE]

    contrast_id <- contrast_row$contrast_id
    write_tsv(result_table, file.path(all_gene_dir, paste0(contrast_id, ".tsv")))
    write_tsv(sig_zero, file.path(sig_zero_dir, paste0(contrast_id, ".tsv")))
    write_tsv(sig_threshold, file.path(sig_threshold_dir, paste0(contrast_id, ".tsv")))

    plot_title <- paste0(population, ": ", contrast_row$contrast_label, " (", design_label, ")")
    save_volcano_plot(result_table, plot_title, file.path(volcano_adjusted_dir, paste0(contrast_id, ".png")), padj_cutoff, lfc_threshold, label_top_n, p_mode = "adjusted")
    save_volcano_plot(result_table, plot_title, file.path(volcano_raw_dir, paste0(contrast_id, ".png")), padj_cutoff, lfc_threshold, label_top_n, p_mode = "raw")
    save_ma_plot(result_table, plot_title, file.path(ma_dir, paste0(contrast_id, ".png")), padj_cutoff, lfc_threshold)

    summaries[[contrast_id]] <- summarise_result_table(
      result_table = result_table,
      population = population,
      contrast_row = contrast_row,
      design_id = design_id,
      design_label = design_label,
      model_formula = model_formula,
      padj_cutoff = padj_cutoff,
      lfc_threshold = lfc_threshold
    )
  }

  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  write_tsv(summary, file.path(result_root, "contrast_summary.tsv"))
  summary
}

iqr_outlier_flags <- function(x, low = TRUE, high = TRUE, multiplier = 3) {
  x <- as.numeric(x)
  q1 <- stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE)
  q3 <- stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)
  iqr <- q3 - q1
  lower <- q1 - multiplier * iqr
  upper <- q3 + multiplier * iqr
  flags <- rep(FALSE, length(x))
  if (low) flags <- flags | (!is.na(x) & x < lower)
  if (high) flags <- flags | (!is.na(x) & x > upper)
  list(flags = flags, lower = lower, upper = upper)
}

write_empty_note <- function(path, lines) {
  writeLines(lines, con = path)
  invisible(path)
}
