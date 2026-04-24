#!/usr/bin/env Rscript

# This script implements two linked tasks for one sorted nuclei population:
#
# 1. Formal technical covariate selection
#    We screen technical variables using:
#    - pairwise collinearity checks
#    - blind PCA on the VST matrix
#    - weighted covariate-to-PC association metrics
#    - multivariable variancePartition
#    The goal is to return an explicit, machine-readable set of technical
#    covariates to carry forward into DEA.
#
# 2. DESeq2 design scouting
#    We then fit a small number of interpretable DESeq2 formulas, including a
#    recommended formula built from the selected covariates, and write the usual
#    PCA / distance / size-factor / dispersion diagnostics for each.
#
# The script assumes the unified handoff  object from 01_build_salmon_gene_inputs.R
# and the EXP383 metadata structure.

suppressPackageStartupMessages({
  required_packages <- c("DESeq2", "ggplot2", "pheatmap", "variancePartition")
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

# Load the shared plotting helpers from this repo so that DEA diagnostics and
# global QC figures use the same population colours and base theme.
get_script_dir <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 0) {
    return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  }

  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE))
}

source(file.path(get_script_dir(), "plot_style.R"))

# -----------------------------
# Command-line argument parsing
# -----------------------------
# This script is intentionally prescriptive, so the command-line arguments are kept
# to a minimum, hard-coding the params below. 
# The only required argument is the population to analyse.
parse_cli_args <- function() {
  defaults <- list(
    input_rds = "results/dea/01_build_salmon_gene_inputs/exp383_salmon_gene_input.rds",
    population = NA_character_,
    design_tsv = "config/dea_designs_initial.tsv",
    output_root = "results/dea/02_design_scout",
    min_count = "10",
    min_samples = "auto",
    pc_variance_threshold = "0.85",
    pc_padj_cutoff = "0.05",
    weighted_pc_cutoff = "0.025",
    varpart_q3_cutoff = "0.03",
    varpart_max_cutoff = "0.75"
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

  numeric_args <- c(
    "min_count",
    "pc_variance_threshold",
    "pc_padj_cutoff",
    "weighted_pc_cutoff",
    "varpart_q3_cutoff",
    "varpart_max_cutoff"
  )
  for (arg_name in numeric_args) {
    parsed[[arg_name]] <- as.numeric(parsed[[arg_name]])
  }

  if (!identical(parsed$min_samples, "auto")) {
    parsed$min_samples <- as.integer(parsed$min_samples)
  }

  parsed
}

# -----------------
# Small utilities
# -----------------
# These helpers are intentionally simple and local so the script remains easy to
# read as a single analysis stage.
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
# defensive preparation here so this script is explicit about what it expects to
# model and plot.
prepare_coldata_for_modeling <- function(coldata) {
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
  scale_map <- c(
    "sample_mass_mg" = "sample_mass_mg_z",
    "incubation_time_hrs" = "incubation_time_hrs_z",
    "date_nuc_prep_days" = "date_nuc_prep_days_z",
    "raw_total_sequences_millions" = "raw_total_sequences_millions_z",
    "raw_percent_duplicates" = "raw_percent_duplicates_z",
    "raw_percent_gc" = "raw_percent_gc_z",
    "fastp_percent_duplication" = "fastp_percent_duplication_z",
    "fastp_q30_rate_after_filtering" = "fastp_q30_rate_after_filtering_z",
    "fastp_passed_filter_reads_millions" = "fastp_passed_filter_reads_millions_z",
    "fastp_percent_gc_after_filtering" = "fastp_percent_gc_after_filtering_z",
    "fastp_percent_surviving" = "fastp_percent_surviving_z",
    "fastp_percent_adapter" = "fastp_percent_adapter_z",
    "rrna_alignment_percent" = "rrna_alignment_percent_z",
    "filtered_percent_duplicates" = "filtered_percent_duplicates_z",
    "filtered_percent_gc" = "filtered_percent_gc_z",
    "filtered_total_sequences_millions" = "filtered_total_sequences_millions_z",
    "salmon_percent_mapped" = "salmon_percent_mapped_z",
    "salmon_num_mapped_millions" = "salmon_num_mapped_millions_z"
  )

  for (source_col in names(scale_map)) {
    target_col <- scale_map[[source_col]]
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
# These helpers are kept in-script because they are tightly coupled to the files
# written by this analysis stage.
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

  exp383_open_png_device(output_file, width = 9, height = 8)
  pheatmap::pheatmap(
    sample_dist_matrix,
    annotation_col = annotation_df,
    annotation_row = annotation_df,
    annotation_colors = exp383_population_annotation_colors(),
    main = "Sample-to-sample distances"
  )
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

save_matrix_heatmap <- function(matrix_data, output_file, title_text, palette_values, breaks = NA) {
  plot_matrix <- matrix_data
  plot_matrix[is.na(plot_matrix)] <- 0

  exp383_open_png_device(output_file, width = 10, height = 8)
  if (all(is.na(breaks))) {
    pheatmap::pheatmap(
      plot_matrix,
      main = title_text,
      color = palette_values
    )
  } else {
    pheatmap::pheatmap(
      plot_matrix,
      main = title_text,
      color = palette_values,
      breaks = breaks
    )
  }
  grDevices::dev.off()
}

# Build a simple matrix from long-format data for pheatmap. This keeps the
# heatmap-writing code readable in the main analysis block.
build_heatmap_matrix <- function(long_table, row_column, column_column, value_column, fill_value = NA_real_) {
  row_values <- unique(long_table[[row_column]])
  column_values <- unique(long_table[[column_column]])

  out <- matrix(fill_value, nrow = length(row_values), ncol = length(column_values))
  rownames(out) <- row_values
  colnames(out) <- column_values

  for (i in seq_len(nrow(long_table))) {
    out[long_table[[row_column]][[i]], long_table[[column_column]][[i]]] <- long_table[[value_column]][[i]]
  }

  out
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
  registry <- data.frame(
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
    stringsAsFactors = FALSE
  )

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
  vp_data <- as.data.frame(coldata[, formal_columns, drop = FALSE])
  assay_matrix <- SummarizedExperiment::assay(vsd)

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
apply_formal_selection_rule <- function(registry, pc_summary, varpart_summary, args) {
  selection_metrics <- merge(registry, pc_summary, by = c("covariate", "analysis_column"), all.x = TRUE)
  selection_metrics <- merge(selection_metrics, varpart_summary, by = c("covariate", "analysis_column"), all.x = TRUE)

  selection_metrics$pc_rule_pass <- FALSE
  selection_metrics$varpart_q3_pass <- !is.na(selection_metrics$varpart_q3_fraction) &
    selection_metrics$varpart_q3_fraction >= args$varpart_q3_cutoff
  selection_metrics$varpart_max_pass <- !is.na(selection_metrics$varpart_max_fraction) &
    selection_metrics$varpart_max_fraction >= args$varpart_max_cutoff

  numeric_rows <- selection_metrics$data_type == "numeric"
  categorical_rows <- selection_metrics$data_type == "categorical"

  selection_metrics$pc_rule_pass[numeric_rows] <- (
    selection_metrics$n_selected_pcs_sig[numeric_rows] >= 1 &
      !is.na(selection_metrics$max_weighted_pc_score[numeric_rows]) &
      selection_metrics$max_weighted_pc_score[numeric_rows] >= args$weighted_pc_cutoff
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

# -----------------
# Main script logic
# -----------------
args <- parse_cli_args()
input_rds <- normalizePath(args$input_rds, winslash = "/", mustWork = FALSE)
design_tsv <- normalizePath(args$design_tsv, winslash = "/", mustWork = FALSE)
stop_if_missing(input_rds, "Prepared input RDS")
stop_if_missing(design_tsv, "Design TSV")

analysis_input <- readRDS(input_rds)
required_objects <- c("txi_gene", "sample_metadata", "gene_annotation")
missing_objects <- setdiff(required_objects, names(analysis_input))
if (length(missing_objects) > 0) {
  stop(
    "Input RDS does not look like the output of 01_build_salmon_gene_inputs.R. ",
    "Missing objects: ", paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

txi_gene <- analysis_input$txi_gene
gene_annotation <- analysis_input$gene_annotation
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

if (identical(args$min_samples, "auto")) {
  group_sizes <- table(coldata_pop$group_assignment)
  min_samples <- as.integer(min(group_sizes))
} else {
  min_samples <- args$min_samples
}

keep_genes <- compute_shared_filter(
  count_matrix = txi_pop$counts,
  grouping = coldata_pop$group_assignment,
  min_count = args$min_count,
  min_samples = min_samples
)

txi_pop_filtered <- subset_txi(
  txi = txi_pop,
  keep_rows = keep_genes,
  keep_cols = rep(TRUE, ncol(txi_pop$counts))
)

shared_root <- file.path(args$output_root, args$population, "shared")
dir.create(shared_root, recursive = TRUE, showWarnings = FALSE)

filter_summary <- data.frame(
  population = args$population,
  n_samples = nrow(coldata_pop),
  n_genes_before_filter = nrow(txi_pop$counts),
  n_genes_after_filter = nrow(txi_pop_filtered$counts),
  min_count = args$min_count,
  min_samples = min_samples,
  filter_group_var = "group_assignment",
  stringsAsFactors = FALSE
)

write_tsv(coldata_pop, file.path(shared_root, "population_metadata.tsv"))
write_tsv(filter_summary, file.path(shared_root, "filter_summary.tsv"))
write_tsv(
  data.frame(gene_id = rownames(txi_pop_filtered$counts), stringsAsFactors = FALSE),
  file.path(shared_root, "filtered_gene_ids.tsv")
)

# ------------------------------
# Shared formal selection stage
# ------------------------------
message("Running shared covariate screening for population: ", args$population)

screen_root <- file.path(shared_root, "covariate_screening")
screen_plot_dir <- file.path(screen_root, "plots")
screen_table_dir <- file.path(screen_root, "tables")
dir.create(screen_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(screen_table_dir, recursive = TRUE, showWarnings = FALSE)

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
selected_pcs <- select_pcs_by_cumulative_variance(percent_var, cumulative_threshold = args$pc_variance_threshold)

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
  pc_padj_cutoff = args$pc_padj_cutoff
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
  args = args
)
selected_covariates <- selection_metrics[selection_metrics$selected, , drop = FALSE]
recommended_formula <- build_recommended_design_formula(selection_metrics)

write_tsv(registry, file.path(screen_table_dir, "candidate_covariates.tsv"))
write_tsv(pairwise_associations, file.path(screen_table_dir, "pairwise_covariate_associations.tsv"))
write_tsv(numeric_collinearity_summary, file.path(screen_table_dir, "technical_numeric_collinearity.tsv"))
write_tsv(formal_registry, file.path(screen_table_dir, "retained_covariates.tsv"))
write_tsv(blind_pca_table, file.path(screen_table_dir, "blind_pca_scores.tsv"))
write_tsv(pc_variance_table, file.path(screen_table_dir, "pc_variance.tsv"))
write_tsv(
  subset(pc_variance_table, selected_for_screen),
  file.path(screen_table_dir, "selected_pcs.tsv")
)
write_tsv(pc_covariate_associations, file.path(screen_table_dir, "pc_covariate_associations.tsv"))
write_tsv(pc_covariate_summary, file.path(screen_table_dir, "pc_covariate_summary.tsv"))
write_tsv(
  variance_partition_results$input_summary,
  file.path(screen_root, "variance_partition_input_summary.tsv")
)
write_tsv(variance_partition_summary, file.path(screen_root, "variance_partition_multivariate.tsv"))
write_tsv(selection_metrics, file.path(screen_root, "covariate_selection_metrics.tsv"))
write_tsv(selected_covariates, file.path(screen_root, "selected_covariates.tsv"))
write_tsv(
  data.frame(
    pc_variance_threshold = args$pc_variance_threshold,
    pc_padj_cutoff = args$pc_padj_cutoff,
    weighted_pc_cutoff = args$weighted_pc_cutoff,
    varpart_q3_cutoff = args$varpart_q3_cutoff,
    varpart_max_cutoff = args$varpart_max_cutoff,
    stringsAsFactors = FALSE
  ),
  file.path(screen_root, "covariate_selection_thresholds.tsv")
)
writeLines(recommended_formula, con = file.path(screen_root, "recommended_design_formula.txt"))
writeLines(
  capture.output(variance_partition_results$varpart),
  con = file.path(screen_root, "variance_partition_model.txt")
)

save_pca_scree_plot(
  pc_variance_table = pc_variance_table,
  n_selected_pcs = length(selected_pcs),
  variance_threshold = args$pc_variance_threshold,
  output_file = file.path(screen_plot_dir, "blind_pca_scree.png")
)

save_numeric_correlation_heatmap(
  coldata = coldata_pop,
  registry = formal_registry,
  output_file = file.path(screen_plot_dir, "retained_numeric_covariate_correlations.png")
)

pc_effect_matrix <- build_heatmap_matrix(
  long_table = transform(pc_covariate_associations, heatmap_value = effect_size_abs),
  row_column = "covariate",
  column_column = "pc",
  value_column = "heatmap_value",
  fill_value = NA_real_
)
save_matrix_heatmap(
  matrix_data = pc_effect_matrix,
  output_file = file.path(screen_plot_dir, "pc_covariate_effect_sizes.png"),
  title_text = "Selected-PC covariate association effect sizes",
  palette_values = grDevices::colorRampPalette(c("white", "#b2182b"))(100)
)

pc_pvalue_matrix <- build_heatmap_matrix(
  long_table = transform(pc_covariate_associations, neg_log10_p = -log10(p_value)),
  row_column = "covariate",
  column_column = "pc",
  value_column = "neg_log10_p",
  fill_value = 0
)
save_matrix_heatmap(
  matrix_data = pc_pvalue_matrix,
  output_file = file.path(screen_plot_dir, "pc_covariate_neglog10_pvalues.png"),
  title_text = "Selected-PC covariate significance (-log10 p)",
  palette_values = grDevices::colorRampPalette(c("white", "#2166ac"))(100)
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
  formal_registry$covariate[formal_registry$include_in_formal_model & formal_registry$is_candidate]
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

  design_dir <- file.path(args$output_root, args$population, design_id)
  plot_dir <- file.path(design_dir, "plots")
  table_dir <- file.path(design_dir, "tables")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  message("Scouting design: ", design_id, " -> ", design_formula)

  writeLines(
    c(
      paste0("population: ", args$population),
      paste0("design_id: ", design_id),
      paste0("design_formula: ", design_formula),
      paste0("min_count: ", args$min_count),
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

  normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
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
  write_tsv(as.data.frame(normalized_counts), file.path(table_dir, "normalized_counts.tsv"), row_names = TRUE)
  write_tsv(pca_table, file.path(table_dir, "pca_coordinates.tsv"))
  write_tsv(as.data.frame(SummarizedExperiment::colData(dds)), file.path(table_dir, "sample_metadata_used.tsv"), row_names = TRUE)
  write_tsv(
    subset(gene_annotation, gene_id %in% rownames(dds)),
    file.path(table_dir, "gene_annotation_used.tsv")
  )
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

  variables_to_plot <- unique(c(
    "group_assignment",
    all.vars(design_formula_obj),
    formal_registry$covariate[formal_registry$is_candidate],
    c("sample_mass_mg", "incubation_time_hrs", "date_nuc_prep_days")
  ))

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

  writeLines(capture.output(sessionInfo()), con = file.path(design_dir, "sessionInfo.txt"))
}

writeLines(capture.output(sessionInfo()), con = file.path(screen_root, "sessionInfo.txt"))

message("Finished DESeq2 design scouting for population: ", args$population)
