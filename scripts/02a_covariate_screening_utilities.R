# Utility functions for scripts/02b_covariate_screening_main_script.R.
# Keep executable workflow code in the main script; keep support code here.

# -----------------------------
# Command-line argument parsing
# -----------------------------
# Only two inputs are exposed at the command line:
# - which sorted population to analyse from start to finish
# - which design table to use for the DESeq2 design-comparison stage
parse_cli_args <- function() {
  defaults <- list(
    population = NA_character_,
    design_tsv = "config/dea_designs_initial.tsv"
  )

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %% 2 != 0) {
    stop(
      "Arguments must be provided as --key value pairs. ",
      "Received an odd number of command-line tokens.",
      call. = FALSE
    )
  }

  parsed <- defaults
  if (length(args) > 0) {
    for (i in seq(1, length(args), by = 2)) {
      key <- args[[i]]
      value <- args[[i + 1]]

      if (!startsWith(key, "--")) {
        stop(sprintf("Unexpected argument name: %s", key), call. = FALSE)
      }

      key <- sub("^--", "", key)
      if (!key %in% names(parsed)) {
        stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
      }

      parsed[[key]] <- value
    }
  }

  if (is.na(parsed$population) || parsed$population == "") {
    stop("--population is required, e.g. --population NeuN", call. = FALSE)
  }

  parsed
}

# -----------------
# Small utilities
# -----------------
# These helpers are intentionally simple so the main workflow can read as a
# single analysis stage.
stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
}

write_tsv <- function(x, path, row_names = FALSE) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = row_names,
    col.names = TRUE
  )
}

sanitize_design_id <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (nchar(x) == 0) {
    x <- "design"
  }
  x
}

get_sample_scaling_factors <- function(dds) {
  size_factors <- DESeq2::sizeFactors(dds)
  if (!is.null(size_factors)) {
    return(size_factors)
  }

  normalization_factors <- DESeq2::normalizationFactors(dds)
  if (is.null(normalization_factors)) {
    stop(
      "DESeq2 object has neither sizeFactors nor normalizationFactors.",
      call. = FALSE
    )
  }

  exp(colMeans(log(normalization_factors)))
}

# -----------------------------
# Metadata and tximport helpers
# -----------------------------
# Script 01 already cleans the unified metadata, but we keep a small amount of
# defensive preparation here so the workflow is explicit about what it expects
# to model and plot.
get_covariate_spec <- function() {
  data.frame(
    covariate = c(
      "group_assignment",
      "inoculation_batch",
      "sample_mass_mg",
      "incubation_time_hrs",
      "date_nuc_prep_days",
      "raw_total_sequences_millions",
      "raw_percent_duplicates",
      "raw_percent_gc",
      "fastp_percent_duplication",
      "fastp_q30_rate_after_filtering",
      "fastp_passed_filter_reads_millions",
      "fastp_percent_gc_after_filtering",
      "fastp_percent_surviving",
      "fastp_percent_adapter",
      "rrna_alignment_percent",
      "filtered_percent_duplicates",
      "filtered_percent_gc",
      "filtered_total_sequences_millions",
      "salmon_percent_mapped",
      "salmon_num_mapped_millions"
    ),
    analysis_column = c(
      "group_assignment",
      "inoculation_batch",
      "sample_mass_mg_z",
      "incubation_time_hrs_z",
      "date_nuc_prep_days_z",
      "raw_total_sequences_millions_z",
      "raw_percent_duplicates_z",
      "raw_percent_gc_z",
      "fastp_percent_duplication_z",
      "fastp_q30_rate_after_filtering_z",
      "fastp_passed_filter_reads_millions_z",
      "fastp_percent_gc_after_filtering_z",
      "fastp_percent_surviving_z",
      "fastp_percent_adapter_z",
      "rrna_alignment_percent_z",
      "filtered_percent_duplicates_z",
      "filtered_percent_gc_z",
      "filtered_total_sequences_millions_z",
      "salmon_percent_mapped_z",
      "salmon_num_mapped_millions_z"
    ),
    role = c(
      "protected_biological",
      "technical_metadata",
      "technical_metadata",
      "technical_metadata",
      "technical_metadata",
      rep("technical_qc", 15)
    ),
    scale_source = c(
      NA_character_,
      NA_character_,
      "sample_mass_mg",
      "incubation_time_hrs",
      "date_nuc_prep_days",
      "raw_total_sequences_millions",
      "raw_percent_duplicates",
      "raw_percent_gc",
      "fastp_percent_duplication",
      "fastp_q30_rate_after_filtering",
      "fastp_passed_filter_reads_millions",
      "fastp_percent_gc_after_filtering",
      "fastp_percent_surviving",
      "fastp_percent_adapter",
      "rrna_alignment_percent",
      "filtered_percent_duplicates",
      "filtered_percent_gc",
      "filtered_total_sequences_millions",
      "salmon_percent_mapped",
      "salmon_num_mapped_millions"
    ),
    stringsAsFactors = FALSE
  )
}

prepare_coldata_for_modeling <- function(coldata) {
  # The prepared input object stores metadata as a tibble. For Bioconductor model
  # code we want an ordinary data.frame with stable sample IDs as row names.
  coldata <- as.data.frame(coldata, stringsAsFactors = FALSE)
  if ("sample" %in% colnames(coldata)) {
    rownames(coldata) <- coldata$sample
  }

  is_character <- vapply(coldata, is.character, logical(1))
  coldata[is_character] <- lapply(coldata[is_character], factor)
  coldata
}

scale_numeric_covariate <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))

  keep <- !is.na(x)
  if (!any(keep)) {
    return(out)
  }

  x_keep <- x[keep]
  if (stats::sd(x_keep) == 0) {
    out[keep] <- 0
    return(out)
  }

  out[keep] <- as.numeric(scale(x_keep))
  out
}

ensure_standardized_covariates <- function(coldata) {
  covariate_spec <- get_covariate_spec()
  scale_rows <- !is.na(covariate_spec$scale_source)

  for (i in which(scale_rows)) {
    source_col <- covariate_spec$scale_source[[i]]
    target_col <- covariate_spec$analysis_column[[i]]
    if (source_col %in% colnames(coldata) && !target_col %in% colnames(coldata)) {
      coldata[[target_col]] <- scale_numeric_covariate(coldata[[source_col]])
    }
  }

  coldata
}

subset_txi <- function(txi, keep_rows, keep_cols) {
  out <- txi
  out$counts <- txi$counts[keep_rows, keep_cols, drop = FALSE]
  out$abundance <- txi$abundance[keep_rows, keep_cols, drop = FALSE]
  out$length <- txi$length[keep_rows, keep_cols, drop = FALSE]
  out$countsFromAbundance <- txi$countsFromAbundance
  out
}

compute_shared_filter <- function(count_matrix, grouping, min_count, min_samples) {
  split_indices <- split(seq_len(ncol(count_matrix)), grouping)

  keep <- rep(FALSE, nrow(count_matrix))
  for (group_name in names(split_indices)) {
    group_idx <- split_indices[[group_name]]
    keep_in_group <- rowSums(count_matrix[, group_idx, drop = FALSE] >= min_count) >= min_samples
    keep <- keep | keep_in_group
  }

  keep
}

# -------------------------
# Plotting helper functions
# -------------------------
# These helpers are tightly coupled to the files written by this analysis stage.
compute_pca_table <- function(vsd, coldata) {
  pca <- stats::prcomp(t(SummarizedExperiment::assay(vsd)), center = TRUE, scale. = FALSE)
  percent_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

  score_table <- as.data.frame(pca$x)
  score_table$sample <- rownames(score_table)

  metadata_df <- as.data.frame(coldata)
  metadata_df$sample <- rownames(coldata)
  metadata_df <- metadata_df[match(score_table$sample, metadata_df$sample), , drop = FALSE]

  extra_columns <- setdiff(colnames(metadata_df), colnames(score_table))
  score_table <- cbind(score_table, metadata_df[, extra_columns, drop = FALSE])
  attr(score_table, "percent_var") <- percent_var
  score_table
}

save_pca_plot <- function(pca_table, color_var, title_text, output_file) {
  if (!color_var %in% colnames(pca_table)) {
    return(invisible(NULL))
  }

  percent_var <- attr(pca_table, "percent_var")
  p <- ggplot2::ggplot(
    pca_table,
    ggplot2::aes(x = PC1, y = PC2, color = .data[[color_var]])
  ) +
    ggplot2::geom_point(size = 2.8, alpha = 0.9) +
    exp383_theme(base_size = 12) +
    ggplot2::labs(
      title = title_text,
      x = sprintf("PC1 (%.2f%% variance)", percent_var[[1]]),
      y = sprintf("PC2 (%.2f%% variance)", percent_var[[2]]),
      color = color_var
    )

  if (identical(color_var, "population")) {
    p <- p + exp383_scale_colour_population(name = "population")
  }

  exp383_save_ggplot(output_file, plot = p, width = 7, height = 5.5)
}

save_sample_distance_heatmap <- function(vsd, annotation_df, output_file) {
  sample_dist_matrix <- as.matrix(stats::dist(t(SummarizedExperiment::assay(vsd))))
  sample_names <- colnames(SummarizedExperiment::assay(vsd))

  # pheatmap expects annotation rows to be named with the same sample IDs as the
  # matrix columns. We also only pass custom annotation colours when the
  # annotation actually includes population, because the shared palette helper is
  # population-specific.
  annotation_df <- as.data.frame(annotation_df)
  if (ncol(annotation_df) > 0) {
    annotation_df <- annotation_df[match(sample_names, rownames(annotation_df)), , drop = FALSE]
    rownames(annotation_df) <- sample_names
  }

  annotation_colors <- NULL
  if ("population" %in% colnames(annotation_df)) {
    annotation_colors <- exp383_population_annotation_colors()
  }

  exp383_open_png_device(output_file, width = 9, height = 8)
  if (ncol(annotation_df) == 0) {
    pheatmap::pheatmap(
      sample_dist_matrix,
      main = "Sample-to-sample distances"
    )
  } else {
    pheatmap::pheatmap(
      sample_dist_matrix,
      annotation_col = annotation_df,
      annotation_row = annotation_df,
      annotation_colors = annotation_colors,
      main = "Sample-to-sample distances"
    )
  }
  grDevices::dev.off()
}

save_dispersion_plot <- function(dds, output_file) {
  exp383_open_png_device(output_file, width = 7, height = 5.5)
  DESeq2::plotDispEsts(dds)
  grDevices::dev.off()
}

save_size_factor_plot <- function(size_factor_table, output_file) {
  p <- ggplot2::ggplot(
    size_factor_table,
    ggplot2::aes(x = reorder(sample, size_factor), y = size_factor, fill = group_assignment)
  ) +
    ggplot2::geom_col(show.legend = TRUE) +
    ggplot2::coord_flip() +
    exp383_theme(base_size = 10) +
    ggplot2::labs(
      title = "DESeq2 size factors",
      x = "Sample",
      y = "Size factor"
    )

  exp383_save_ggplot(output_file, plot = p, width = 7, height = 10)
}

save_pca_scree_plot <- function(pc_variance_table, n_selected_pcs, variance_threshold, output_file) {
  threshold_percent <- 100 * variance_threshold
  p <- ggplot2::ggplot(pc_variance_table, ggplot2::aes(x = pc_index, y = variance_percent)) +
    ggplot2::geom_col(fill = "#3b6fb6") +
    ggplot2::geom_line(ggplot2::aes(y = cumulative_variance_percent), color = "#b22222", group = 1) +
    ggplot2::geom_point(ggplot2::aes(y = cumulative_variance_percent), color = "#b22222") +
    ggplot2::geom_vline(xintercept = n_selected_pcs + 0.5, linetype = 2, color = "grey40") +
    exp383_theme(base_size = 12) +
    ggplot2::scale_x_continuous(breaks = pc_variance_table$pc_index) +
    ggplot2::labs(
      title = "Blind PCA scree plot",
      subtitle = sprintf(
        "Selected PCs explain at least %.0f%% cumulative variance (%d PCs retained)",
        threshold_percent,
        n_selected_pcs
      ),
      x = "Principal component",
      y = "Variance explained (%)"
    )

  exp383_save_ggplot(output_file, plot = p, width = 8, height = 5)
}

save_numeric_correlation_heatmap <- function(coldata, registry, output_file) {
  numeric_covariates <- registry$analysis_column[
    registry$include_in_formal_model &
      registry$is_candidate &
      registry$data_type == "numeric"
  ]

  if (length(numeric_covariates) < 2) {
    return(invisible(NULL))
  }

  correlation_matrix <- stats::cor(
    as.data.frame(coldata[, numeric_covariates, drop = FALSE]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )

  exp383_open_png_device(output_file, width = 9, height = 8)
  pheatmap::pheatmap(
    correlation_matrix,
    main = "Spearman correlations among retained technical numeric covariates",
    color = grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
    breaks = seq(-1, 1, length.out = 101)
  )
  grDevices::dev.off()
}

save_metric_barplot <- function(plot_table, metric_column, title_text, y_label, output_file) {
  plot_table <- plot_table[!is.na(plot_table[[metric_column]]), , drop = FALSE]
  if (nrow(plot_table) == 0) {
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(
    plot_table,
    ggplot2::aes(x = reorder(covariate, .data[[metric_column]]), y = .data[[metric_column]], fill = selected)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    exp383_theme(base_size = 11) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#16a829", "FALSE" = "#bdbdbd")) +
    ggplot2::labs(
      title = title_text,
      x = "Covariate",
      y = y_label,
      fill = "Selected"
    )

  exp383_save_ggplot(output_file, plot = p, width = 8, height = 6)
}

save_variance_partition_plot <- function(varpart_object, output_file) {
  exp383_open_png_device(output_file, width = 9, height = 7)
  plot_object <- variancePartition::sortCols(varpart_object) |>
    variancePartition::plotVarPart(label.angle = 60)
  print(plot_object)
  grDevices::dev.off()
}

# ----------------------------
# Covariate helper functions
# ----------------------------
# The registry is the central reference table for the screen. It records which
# covariates are protected biological terms, which are technical candidates, and
# which metadata column should actually be modelled.
build_covariate_registry <- function(coldata) {
  registry <- get_covariate_spec()[, c("covariate", "analysis_column", "role")]

  keep_rows <- registry$covariate %in% colnames(coldata) & registry$analysis_column %in% colnames(coldata)
  registry <- registry[keep_rows, , drop = FALSE]
  registry$is_candidate <- registry$role != "protected_biological"
  registry$data_type <- vapply(
    registry$analysis_column,
    function(column_name) {
      if (is.numeric(coldata[[column_name]]) || is.integer(coldata[[column_name]])) {
        "numeric"
      } else {
        "categorical"
      }
    },
    character(1)
  )
  registry$n_non_missing <- vapply(
    registry$analysis_column,
    function(column_name) sum(!is.na(coldata[[column_name]])),
    integer(1)
  )
  registry$n_unique <- vapply(
    registry$analysis_column,
    function(column_name) length(unique(stats::na.omit(coldata[[column_name]]))),
    integer(1)
  )
  registry$n_missing <- nrow(coldata) - registry$n_non_missing
  registry
}

cramers_v <- function(chisq_statistic, n, nrow_table, ncol_table) {
  if (is.na(chisq_statistic) || n <= 0 || min(nrow_table, ncol_table) <= 1) {
    return(NA_real_)
  }

  sqrt(as.numeric(chisq_statistic) / (n * min(nrow_table - 1, ncol_table - 1)))
}

kruskal_epsilon_squared <- function(kruskal_statistic, n, k_groups) {
  if (is.na(kruskal_statistic) || n <= k_groups || k_groups <= 1) {
    return(NA_real_)
  }

  max(0, (as.numeric(kruskal_statistic) - k_groups + 1) / (n - k_groups))
}

# This one function handles all pairwise association testing used in the screen.
# It keeps the output schema stable regardless of whether the pair is numeric /
# numeric, numeric / categorical, or categorical / categorical.
pairwise_covariate_association <- function(x, y, x_name, y_name, x_type, y_type) {
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]
  n_complete <- length(x)

  empty_result <- data.frame(
    x_name = x_name,
    y_name = y_name,
    x_type = x_type,
    y_type = y_type,
    n_complete = n_complete,
    test_name = NA_character_,
    effect_size_name = NA_character_,
    effect_size = NA_real_,
    effect_size_abs = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  if (n_complete < 3) {
    return(empty_result)
  }

  if (x_type == "numeric" && y_type == "numeric") {
    if (length(unique(x)) < 2 || length(unique(y)) < 2) {
      return(empty_result)
    }

    test <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
    empty_result$test_name <- "spearman"
    empty_result$effect_size_name <- "spearman_rho"
    empty_result$effect_size <- unname(test$estimate)
    empty_result$effect_size_abs <- abs(empty_result$effect_size)
    empty_result$statistic <- unname(test$statistic)
    empty_result$p_value <- test$p.value
    return(empty_result)
  }

  if (x_type != y_type) {
    numeric_values <- if (x_type == "numeric") x else y
    grouping <- if (x_type == "categorical") as.factor(x) else as.factor(y)

    if (length(unique(numeric_values)) < 2 || nlevels(grouping) < 2) {
      return(empty_result)
    }

    kw_data <- data.frame(value = numeric_values, group = grouping)
    test <- stats::kruskal.test(value ~ group, data = kw_data)
    empty_result$test_name <- "kruskal_wallis"
    empty_result$effect_size_name <- "epsilon_squared"
    empty_result$effect_size <- kruskal_epsilon_squared(test$statistic, n = n_complete, k_groups = nlevels(grouping))
    empty_result$effect_size_abs <- abs(empty_result$effect_size)
    empty_result$statistic <- unname(test$statistic)
    empty_result$p_value <- test$p.value
    return(empty_result)
  }

  contingency_table <- table(as.factor(x), as.factor(y))
  if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2) {
    return(empty_result)
  }

  test <- suppressWarnings(stats::chisq.test(contingency_table))
  empty_result$test_name <- "chi_square"
  empty_result$effect_size_name <- "cramers_v"
  empty_result$effect_size <- cramers_v(
    chisq_statistic = test$statistic,
    n = sum(contingency_table),
    nrow_table = nrow(contingency_table),
    ncol_table = ncol(contingency_table)
  )
  empty_result$effect_size_abs <- abs(empty_result$effect_size)
  empty_result$statistic <- unname(test$statistic)
  empty_result$p_value <- test$p.value
  empty_result
}

build_pairwise_covariate_associations <- function(coldata, registry) {
  if (nrow(registry) < 2) {
    return(data.frame())
  }

  pairs <- utils::combn(seq_len(nrow(registry)), 2, simplify = FALSE)
  results <- lapply(pairs, function(pair_idx) {
    x_row <- registry[pair_idx[[1]], , drop = FALSE]
    y_row <- registry[pair_idx[[2]], , drop = FALSE]

    pairwise_covariate_association(
      x = coldata[[x_row$analysis_column]],
      y = coldata[[y_row$analysis_column]],
      x_name = x_row$covariate,
      y_name = y_row$covariate,
      x_type = x_row$data_type,
      y_type = y_row$data_type
    )
  })

  do.call(rbind, results)
}

# Numeric collinearity is the only part of the screen that is pruned
# automatically. We keep the variable with the lower overall correlation burden
# so the formal multivariable model is less likely to be unstable.
build_numeric_collinearity_summary <- function(pairwise_table, registry, rho_threshold = 0.75, p_threshold = 0.05) {
  technical_numeric <- registry$covariate[
    registry$is_candidate &
      registry$data_type == "numeric"
  ]

  summary_table <- data.frame(
    covariate = technical_numeric,
    mean_abs_spearman_to_other_technical_numeric = NA_real_,
    retain_for_formal_model = TRUE,
    drop_reason = NA_character_,
    stringsAsFactors = FALSE
  )

  if (length(technical_numeric) == 0) {
    return(summary_table)
  }

  numeric_pairs <- subset(
    pairwise_table,
    x_name %in% technical_numeric &
      y_name %in% technical_numeric &
      test_name == "spearman"
  )

  summary_table$mean_abs_spearman_to_other_technical_numeric <- vapply(
    technical_numeric,
    function(covariate_name) {
      pair_rows <- numeric_pairs$x_name == covariate_name | numeric_pairs$y_name == covariate_name
      values <- numeric_pairs$effect_size_abs[pair_rows]
      if (length(values) == 0 || all(is.na(values))) {
        return(NA_real_)
      }
      mean(values, na.rm = TRUE)
    },
    numeric(1)
  )

  high_pairs <- subset(
    numeric_pairs,
    !is.na(effect_size_abs) &
      effect_size_abs >= rho_threshold &
      !is.na(p_value) &
      p_value <= p_threshold
  )

  if (nrow(high_pairs) == 0) {
    return(summary_table)
  }

  high_pairs <- high_pairs[order(high_pairs$effect_size_abs, decreasing = TRUE), , drop = FALSE]
  for (i in seq_len(nrow(high_pairs))) {
    x_name <- high_pairs$x_name[[i]]
    y_name <- high_pairs$y_name[[i]]
    x_idx <- match(x_name, summary_table$covariate)
    y_idx <- match(y_name, summary_table$covariate)

    if (!summary_table$retain_for_formal_model[[x_idx]] || !summary_table$retain_for_formal_model[[y_idx]]) {
      next
    }

    x_score <- summary_table$mean_abs_spearman_to_other_technical_numeric[[x_idx]]
    y_score <- summary_table$mean_abs_spearman_to_other_technical_numeric[[y_idx]]

    drop_idx <- if (is.na(x_score) && !is.na(y_score)) {
      x_idx
    } else if (!is.na(x_score) && is.na(y_score)) {
      y_idx
    } else if (isTRUE(x_score > y_score)) {
      x_idx
    } else if (isTRUE(y_score > x_score)) {
      y_idx
    } else {
      max(x_idx, y_idx)
    }

    keep_name <- if (drop_idx == x_idx) y_name else x_name
    summary_table$retain_for_formal_model[[drop_idx]] <- FALSE
    summary_table$drop_reason[[drop_idx]] <- sprintf(
      "Dropped after high collinearity with %s (|rho|=%.3f). %s had the lower overall mean absolute correlation to other technical numeric covariates.",
      keep_name,
      high_pairs$effect_size_abs[[i]],
      keep_name
    )
  }

  summary_table
}

finalize_covariate_registry <- function(registry, numeric_collinearity_summary, n_samples) {
  registry$include_in_formal_model <- FALSE
  registry$retention_reason <- "Not assessed."

  for (i in seq_len(nrow(registry))) {
    if (registry$role[[i]] == "protected_biological") {
      registry$include_in_formal_model[[i]] <- TRUE
      registry$retention_reason[[i]] <- "Protected biological term: always included in the formal model."
      next
    }

    if (registry$n_non_missing[[i]] < n_samples) {
      registry$retention_reason[[i]] <- "Dropped because the covariate contains missing values."
      next
    }

    if (registry$n_unique[[i]] < 2) {
      registry$retention_reason[[i]] <- "Dropped because the covariate has no usable variation."
      next
    }

    registry$include_in_formal_model[[i]] <- TRUE
    registry$retention_reason[[i]] <- "Retained for formal selection."
  }

  if (nrow(numeric_collinearity_summary) > 0) {
    idx <- match(numeric_collinearity_summary$covariate, registry$covariate)
    registry$include_in_formal_model[idx] <- numeric_collinearity_summary$retain_for_formal_model
    registry$retention_reason[idx] <- ifelse(
      is.na(numeric_collinearity_summary$drop_reason),
      "Retained after numeric collinearity screen.",
      numeric_collinearity_summary$drop_reason
    )
  }

  registry
}

# -----------------------------
# PCA screening helper methods
# -----------------------------
# PCA is used here as a sample-level screen. We keep PCs until a cumulative
# explained-variance threshold is reached rather than picking an arbitrary fixed
# count.
select_pcs_by_cumulative_variance <- function(percent_var, cumulative_threshold = 0.85) {
  threshold_percent <- 100 * cumulative_threshold
  cumulative_percent <- cumsum(percent_var)
  selected_n <- which(cumulative_percent >= threshold_percent)[1]

  if (is.na(selected_n)) {
    selected_n <- length(percent_var)
  }

  selected_n <- max(2, selected_n)
  paste0("PC", seq_len(selected_n))
}

build_pc_covariate_associations <- function(pc_table, registry, selected_pcs) {
  percent_var <- attr(pc_table, "percent_var")
  formal_covariates <- registry$covariate[registry$include_in_formal_model]

  results <- list()
  result_index <- 1L
  for (pc_name in selected_pcs) {
    pc_index <- as.integer(sub("^PC", "", pc_name))
    pc_variance_fraction <- percent_var[[pc_index]] / 100

    for (covariate_name in formal_covariates) {
      covariate_row <- registry[match(covariate_name, registry$covariate), , drop = FALSE]

      association <- pairwise_covariate_association(
        x = pc_table[[pc_name]],
        y = pc_table[[covariate_row$analysis_column]],
        x_name = pc_name,
        y_name = covariate_name,
        x_type = "numeric",
        y_type = covariate_row$data_type
      )

      association$pc <- pc_name
      association$covariate <- covariate_name
      association$analysis_column <- covariate_row$analysis_column
      association$covariate_role <- covariate_row$role
      association$covariate_type <- covariate_row$data_type
      association$pc_variance_fraction <- pc_variance_fraction
      association$weighted_pc_score <- if (covariate_row$data_type == "numeric") {
        pc_variance_fraction * (association$effect_size_abs ^ 2)
      } else {
        NA_real_
      }

      results[[result_index]] <- association
      result_index <- result_index + 1L
    }
  }

  association_table <- do.call(rbind, results)
  association_table$padj <- ave(
    association_table$p_value,
    association_table$pc,
    FUN = function(x) stats::p.adjust(x, method = "fdr")
  )
  association_table
}

summarize_pc_covariate_associations <- function(pc_association_table, registry, pc_padj_cutoff) {
  formal_covariates <- registry$covariate[registry$include_in_formal_model]

  summary_rows <- lapply(formal_covariates, function(covariate_name) {
    covariate_row <- registry[match(covariate_name, registry$covariate), , drop = FALSE]
    covariate_hits <- pc_association_table[pc_association_table$covariate == covariate_name, , drop = FALSE]

    if (nrow(covariate_hits) == 0) {
      return(data.frame(
        covariate = covariate_name,
        analysis_column = covariate_row$analysis_column,
        covariate_role = covariate_row$role,
        covariate_type = covariate_row$data_type,
        n_selected_pcs = 0,
        n_selected_pcs_sig = 0,
        min_selected_pc_padj = NA_real_,
        max_effect_size_abs = NA_real_,
        max_weighted_pc_score = NA_real_,
        strongest_pc = NA_character_,
        stringsAsFactors = FALSE
      ))
    }

    weighted_values <- covariate_hits$weighted_pc_score
    strongest_pc <- if (all(is.na(weighted_values))) {
      covariate_hits$pc[[which.max(covariate_hits$effect_size_abs)]]
    } else {
      covariate_hits$pc[[which.max(weighted_values)]]
    }

    data.frame(
      covariate = covariate_name,
      analysis_column = covariate_row$analysis_column,
      covariate_role = covariate_row$role,
      covariate_type = covariate_row$data_type,
      n_selected_pcs = nrow(covariate_hits),
      n_selected_pcs_sig = sum(covariate_hits$padj <= pc_padj_cutoff, na.rm = TRUE),
      min_selected_pc_padj = min(covariate_hits$padj, na.rm = TRUE),
      max_effect_size_abs = max(covariate_hits$effect_size_abs, na.rm = TRUE),
      max_weighted_pc_score = if (all(is.na(weighted_values))) NA_real_ else max(weighted_values, na.rm = TRUE),
      strongest_pc = strongest_pc,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summary_rows)
}

# ----------------------------------
# variancePartition helper functions
# ----------------------------------
# The formal gene-level screen is multivariable rather than univariate. This is
# the key shift from the earlier exploratory version of the script.
run_multivariate_variance_partition <- function(vsd, registry, coldata) {
  formal_registry <- registry[registry$include_in_formal_model, , drop = FALSE]
  formal_columns <- formal_registry$analysis_column

  vp_formula <- stats::reformulate(formal_columns)
  assay_matrix <- SummarizedExperiment::assay(vsd)
  sample_names <- colnames(assay_matrix)
  vp_data <- as.data.frame(coldata[, formal_columns, drop = FALSE])
  vp_data <- vp_data[match(sample_names, rownames(vp_data)), , drop = FALSE]
  rownames(vp_data) <- sample_names

  input_summary <- data.frame(
    genes_used_for_variance_partition = nrow(assay_matrix),
    n_covariates_in_formal_model = length(formal_columns),
    formula = paste(deparse(vp_formula), collapse = " "),
    stringsAsFactors = FALSE
  )

  varpart_object <- variancePartition::fitExtractVarPartModel(
    exprObj = assay_matrix,
    formula = vp_formula,
    data = vp_data
  )

  list(
    varpart = varpart_object,
    input_summary = input_summary
  )
}

summarize_variance_partition <- function(varpart_object, registry) {
  formal_registry <- registry[registry$include_in_formal_model, , drop = FALSE]

  summary_rows <- lapply(seq_len(nrow(formal_registry)), function(i) {
    covariate_name <- formal_registry$covariate[[i]]
    analysis_column <- formal_registry$analysis_column[[i]]

    if (!analysis_column %in% colnames(varpart_object)) {
      return(data.frame(
        covariate = covariate_name,
        analysis_column = analysis_column,
        varpart_mean_fraction = NA_real_,
        varpart_median_fraction = NA_real_,
        varpart_q3_fraction = NA_real_,
        varpart_max_fraction = NA_real_,
        varpart_mean_percent = NA_real_,
        varpart_median_percent = NA_real_,
        varpart_q3_percent = NA_real_,
        varpart_max_percent = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    values <- varpart_object[[analysis_column]]
    data.frame(
      covariate = covariate_name,
      analysis_column = analysis_column,
      varpart_mean_fraction = mean(values, na.rm = TRUE),
      varpart_median_fraction = stats::median(values, na.rm = TRUE),
      varpart_q3_fraction = as.numeric(stats::quantile(values, 0.75, na.rm = TRUE)),
      varpart_max_fraction = max(values, na.rm = TRUE),
      varpart_mean_percent = 100 * mean(values, na.rm = TRUE),
      varpart_median_percent = 100 * stats::median(values, na.rm = TRUE),
      varpart_q3_percent = 100 * as.numeric(stats::quantile(values, 0.75, na.rm = TRUE)),
      varpart_max_percent = 100 * max(values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summary_rows)
}

# --------------------------------------
# Formal covariate-selection decision
# --------------------------------------
# Continuous and categorical covariates are handled slightly differently, which
# matches the logic in the Brenton reference scripts:
# - continuous terms must align with important PCs and also show broad gene-level
#   variance via variancePartition
# - categorical terms must show at least one significant PC association and the
#   same broad gene-level variance
# - any term can still be rescued if it has an extreme gene-level effect
apply_formal_selection_rule <- function(registry, pc_summary, varpart_summary, settings) {
  selection_metrics <- merge(registry, pc_summary, by = c("covariate", "analysis_column"), all.x = TRUE)
  selection_metrics <- merge(selection_metrics, varpart_summary, by = c("covariate", "analysis_column"), all.x = TRUE)

  selection_metrics$pc_rule_pass <- FALSE
  selection_metrics$varpart_q3_pass <- !is.na(selection_metrics$varpart_q3_fraction) &
    selection_metrics$varpart_q3_fraction >= settings$varpart_q3_cutoff
  selection_metrics$varpart_max_pass <- !is.na(selection_metrics$varpart_max_fraction) &
    selection_metrics$varpart_max_fraction >= settings$varpart_max_cutoff

  numeric_rows <- selection_metrics$data_type == "numeric"
  categorical_rows <- selection_metrics$data_type == "categorical"

  selection_metrics$pc_rule_pass[numeric_rows] <- (
    selection_metrics$n_selected_pcs_sig[numeric_rows] >= 1 &
      !is.na(selection_metrics$max_weighted_pc_score[numeric_rows]) &
      selection_metrics$max_weighted_pc_score[numeric_rows] >= settings$weighted_pc_cutoff
  )

  selection_metrics$pc_rule_pass[categorical_rows] <- (
    selection_metrics$n_selected_pcs_sig[categorical_rows] >= 1
  )

  selection_metrics$selected <- (
    selection_metrics$is_candidate &
      (
        (selection_metrics$pc_rule_pass & selection_metrics$varpart_q3_pass) |
          selection_metrics$varpart_max_pass
      )
  )

  selection_metrics$selection_reason <- ifelse(
    !selection_metrics$is_candidate,
    "Protected biological term: included in the formal model but not part of technical covariate selection.",
    ifelse(
      selection_metrics$selected & selection_metrics$varpart_max_pass,
      "Selected because the maximum gene-level variance explained exceeded the rescue cutoff.",
      ifelse(
        selection_metrics$selected,
        "Selected because the PC-association rule and variancePartition Q3 cutoff were both met.",
        "Not selected by the formal covariate rule."
      )
    )
  )

  selection_metrics
}

build_recommended_design_formula <- function(selection_metrics) {
  selected_model_columns <- selection_metrics$analysis_column[selection_metrics$selected]
  design_terms <- unique(c("group_assignment", selected_model_columns))
  paste("~", paste(design_terms, collapse = " + "))
}

load_design_table <- function(design_tsv) {
  stop_if_missing(design_tsv, "Design TSV")

  design_table <- read.delim(design_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  required_columns <- c("design_id", "design_formula")
  missing_columns <- setdiff(required_columns, colnames(design_table))
  if (length(missing_columns) > 0) {
    stop(
      sprintf("Design TSV must contain columns: %s", paste(required_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  design_table$design_id <- vapply(design_table$design_id, sanitize_design_id, character(1))
  design_table
}

append_recommended_design <- function(design_table, recommended_formula) {
  recommended_id <- "selected_covariates"
  if (any(design_table$design_formula == recommended_formula)) {
    return(design_table)
  }

  rbind(
    data.frame(
      design_id = recommended_id,
      design_formula = recommended_formula,
      stringsAsFactors = FALSE
    ),
    design_table
  )
}
