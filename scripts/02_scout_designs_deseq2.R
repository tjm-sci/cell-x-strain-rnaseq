#!/usr/bin/env Rscript

# This script is the design-scouting stage for DEA.
#
# It sits between the stable tximport handoff built by
# 01_build_salmon_gene_inputs.R and any final inferential DEA script.
#
# The workflow has two layers:
# 1. A shared, population-level covariate screen inspired by the paper the user
#    cited. We first reduce obviously redundant technical covariates, then ask
#    which remaining variables line up with the dominant PCs of the expression
#    matrix.
# 2. A design-specific DESeq2 QC pass. Once we have a shortlist of plausible
#    covariates, we still want to inspect how candidate formulas behave when we
#    actually fit DESeq2 models.
#
# Output layout:
#   results/dea/02_design_scout/<population>/shared
#   results/dea/02_design_scout/<population>/shared/covariate_screening
#   results/dea/02_design_scout/<population>/<design_id>/plots
#   results/dea/02_design_scout/<population>/<design_id>/tables

suppressPackageStartupMessages({
  required_packages <- c("DESeq2", "ggplot2", "pheatmap")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      ". Install them before running this script.",
      call. = FALSE
    )
  }
})

# Load the shared plot helpers from the same scripts directory so all DEA/QC
# plots use the same palette and base theme as the FANS analysis repo.
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
parse_cli_args <- function() {
  defaults <- list(
    input_rds = "results/dea/01_build_salmon_gene_inputs/exp383_salmon_gene_input.rds",
    population = NA_character_,
    design_tsv = NA_character_,
    design_id = NA_character_,
    design_formula = NA_character_,
    output_root = "results/dea/02_design_scout",
    min_count = "10",
    min_samples = "auto",
    filter_group_var = "group_assignment"
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

  using_design_tsv <- !is.na(parsed$design_tsv) && parsed$design_tsv != ""
  using_single_design <- !is.na(parsed$design_id) && !is.na(parsed$design_formula) &&
    parsed$design_id != "" && parsed$design_formula != ""

  if (!using_design_tsv && !using_single_design) {
    stop(
      "Provide either --design_tsv <file> or both --design_id and --design_formula.",
      call. = FALSE
    )
  }

  parsed$min_count <- as.numeric(parsed$min_count)
  if (!identical(parsed$min_samples, "auto")) {
    parsed$min_samples <- as.integer(parsed$min_samples)
  }

  parsed
}

# -----------------
# Utility functions
# -----------------
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

load_design_table <- function(args) {
  if (!is.na(args$design_tsv) && args$design_tsv != "") {
    stop_if_missing(args$design_tsv, "Design TSV")
    design_table <- read.delim(args$design_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    required_columns <- c("design_id", "design_formula")
    missing_columns <- setdiff(required_columns, colnames(design_table))
    if (length(missing_columns) > 0) {
      stop(
        sprintf(
          "Design TSV must contain columns: %s",
          paste(required_columns, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  } else {
    design_table <- data.frame(
      design_id = args$design_id,
      design_formula = args$design_formula,
      stringsAsFactors = FALSE
    )
  }

  design_table$design_id <- vapply(design_table$design_id, sanitize_design_id, character(1))
  design_table
}

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

# Script 01 now creates the standardised numeric columns, but we also make the
# design scout robust to older handoff objects that may not yet contain them.
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

compute_shared_filter <- function(count_matrix, coldata, min_count, min_samples, filter_group_var) {
  if (!filter_group_var %in% colnames(coldata)) {
    stop(sprintf("filter_group_var '%s' is not present in metadata.", filter_group_var), call. = FALSE)
  }

  grouping <- coldata[[filter_group_var]]
  split_indices <- split(seq_len(ncol(count_matrix)), grouping)

  keep <- rep(FALSE, nrow(count_matrix))
  for (group_name in names(split_indices)) {
    group_idx <- split_indices[[group_name]]
    group_hits <- rowSums(count_matrix[, group_idx, drop = FALSE] >= min_count) >= min_samples
    keep <- keep | group_hits
  }

  keep
}

subset_txi <- function(txi, keep_rows, keep_cols) {
  out <- txi
  out$counts <- txi$counts[keep_rows, keep_cols, drop = FALSE]
  out$abundance <- txi$abundance[keep_rows, keep_cols, drop = FALSE]
  out$length <- txi$length[keep_rows, keep_cols, drop = FALSE]
  out$countsFromAbundance <- txi$countsFromAbundance
  out
}

build_pca_table <- function(vsd, coldata) {
  pca <- DESeq2::plotPCA(vsd, intgroup = "group_assignment", returnData = TRUE)
  pca$sample <- rownames(pca)

  metadata_df <- as.data.frame(coldata)
  metadata_df$sample <- rownames(coldata)
  metadata_df <- metadata_df[match(pca$sample, metadata_df$sample), , drop = FALSE]

  extra_columns <- setdiff(colnames(metadata_df), colnames(pca))
  cbind(pca, metadata_df[, extra_columns, drop = FALSE])
}

save_pca_plot <- function(pca_table, color_var, title_text, output_file) {
  if (!color_var %in% colnames(pca_table)) {
    return(invisible(NULL))
  }

  plot_data <- pca_table
  percent_var <- attr(pca_table, "percent_var")
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC1, y = PC2, color = .data[[color_var]])) +
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

  ggplot2::ggsave(output_file, plot = p, width = 7, height = 5.5)
}

save_sample_distance_heatmap <- function(vsd, annotation_df, output_file) {
  sample_dist <- dist(t(SummarizedExperiment::assay(vsd)))
  sample_dist_matrix <- as.matrix(sample_dist)

  grDevices::pdf(output_file, width = 9, height = 8)
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
  grDevices::pdf(output_file, width = 7, height = 5.5)
  DESeq2::plotDispEsts(dds)
  grDevices::dev.off()
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

  ggplot2::ggsave(output_file, plot = p, width = 7, height = 10)
}

count_non_missing_unique <- function(x) {
  unique(stats::na.omit(x))
}

build_covariate_registry <- function(coldata) {
  registry <- data.frame(
    covariate = c(
      "group_assignment",
      "inoculum",
      "dpi_factor",
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
    role = c(
      "biological",
      "biological",
      "biological",
      "technical_metadata",
      "technical_metadata",
      "technical_metadata",
      "technical_metadata",
      rep("technical_qc", 15)
    ),
    model_column = c(
      "group_assignment",
      "inoculum",
      "dpi_factor",
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
    stringsAsFactors = FALSE
  )

  registry <- registry[registry$covariate %in% colnames(coldata), , drop = FALSE]
  registry$data_type <- vapply(
    registry$covariate,
    function(covariate_name) {
      if (is.numeric(coldata[[covariate_name]]) || is.integer(coldata[[covariate_name]])) {
        "numeric"
      } else {
        "categorical"
      }
    },
    character(1)
  )
  registry$n_non_missing <- vapply(
    registry$covariate,
    function(covariate_name) sum(!is.na(coldata[[covariate_name]])),
    integer(1)
  )
  registry$n_unique <- vapply(
    registry$covariate,
    function(covariate_name) length(count_non_missing_unique(coldata[[covariate_name]])),
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

  pairs <- utils::combn(registry$covariate, 2, simplify = FALSE)
  results <- lapply(pairs, function(pair_names) {
    x_name <- pair_names[[1]]
    y_name <- pair_names[[2]]
    x_type <- registry$data_type[match(x_name, registry$covariate)]
    y_type <- registry$data_type[match(y_name, registry$covariate)]

    pairwise_covariate_association(
      x = coldata[[x_name]],
      y = coldata[[y_name]],
      x_name = x_name,
      y_name = y_name,
      x_type = x_type,
      y_type = y_type
    )
  })

  do.call(rbind, results)
}

build_numeric_collinearity_summary <- function(pairwise_table, registry, rho_threshold = 0.75, p_threshold = 0.05) {
  technical_numeric <- registry$covariate[registry$role != "biological" & registry$data_type == "numeric"]
  if (length(technical_numeric) == 0) {
    return(data.frame())
  }

  numeric_pairs <- subset(
    pairwise_table,
    x_name %in% technical_numeric &
      y_name %in% technical_numeric &
      test_name == "spearman"
  )

  mean_abs_cor <- vapply(technical_numeric, function(covariate_name) {
    pair_rows <- numeric_pairs$x_name == covariate_name | numeric_pairs$y_name == covariate_name
    values <- numeric_pairs$effect_size_abs[pair_rows]
    if (length(values) == 0 || all(is.na(values))) {
      return(NA_real_)
    }
    mean(values, na.rm = TRUE)
  }, numeric(1))

  summary_table <- data.frame(
    covariate = technical_numeric,
    mean_abs_spearman_to_other_technical_numeric = mean_abs_cor,
    retain_for_pc_screening = TRUE,
    drop_reason = NA_character_,
    stringsAsFactors = FALSE
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

    if (!summary_table$retain_for_pc_screening[[x_idx]] || !summary_table$retain_for_pc_screening[[y_idx]]) {
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
      # Deterministic tie-break so reruns are stable.
      max(x_idx, y_idx)
    }

    keep_name <- if (drop_idx == x_idx) y_name else x_name
    summary_table$retain_for_pc_screening[[drop_idx]] <- FALSE
    summary_table$drop_reason[[drop_idx]] <- sprintf(
      "Dropped after high collinearity with %s (|rho|=%.3f). %s had the lower overall mean absolute correlation to other technical numeric covariates.",
      keep_name,
      high_pairs$effect_size_abs[[i]],
      keep_name
    )
  }

  summary_table[order(summary_table$retain_for_pc_screening, summary_table$covariate, decreasing = TRUE), , drop = FALSE]
}

build_retained_covariate_table <- function(registry, numeric_collinearity_summary) {
  retained <- registry
  retained$retain_for_pc_screening <- TRUE
  retained$retention_reason <- ifelse(
    retained$role == "biological",
    "Always retained for shared screening as a biological variable.",
    "Retained."
  )

  if (nrow(numeric_collinearity_summary) > 0) {
    idx <- match(numeric_collinearity_summary$covariate, retained$covariate)
    retained$retain_for_pc_screening[idx] <- numeric_collinearity_summary$retain_for_pc_screening
    drop_reasons <- numeric_collinearity_summary$drop_reason
    retained$retention_reason[idx] <- ifelse(
      is.na(drop_reasons),
      "Retained after technical numeric collinearity screen.",
      drop_reasons
    )
  }

  retained
}

build_pc_score_table <- function(vsd, coldata, n_pcs = 10) {
  expression_matrix <- t(SummarizedExperiment::assay(vsd))
  pca <- stats::prcomp(expression_matrix, center = TRUE, scale. = FALSE)
  n_keep <- min(n_pcs, ncol(pca$x))

  pc_scores <- as.data.frame(pca$x[, seq_len(n_keep), drop = FALSE])
  colnames(pc_scores) <- paste0("PC", seq_len(n_keep))
  pc_scores$sample <- rownames(pc_scores)

  metadata_df <- as.data.frame(coldata)
  metadata_df$sample <- rownames(coldata)
  metadata_df <- metadata_df[match(pc_scores$sample, metadata_df$sample), , drop = FALSE]

  extra_columns <- setdiff(colnames(metadata_df), colnames(pc_scores))
  pc_scores <- cbind(pc_scores, metadata_df[, extra_columns, drop = FALSE])

  percent_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  attr(pc_scores, "percent_var") <- percent_var[seq_len(n_keep)]
  pc_scores
}

build_pc_covariate_associations <- function(pc_table, retained_registry) {
  pc_columns <- grep("^PC[0-9]+$", colnames(pc_table), value = TRUE)
  retained_covariates <- retained_registry$covariate[retained_registry$retain_for_pc_screening]

  results <- list()
  result_index <- 1L
  for (pc_name in pc_columns) {
    for (covariate_name in retained_covariates) {
      covariate_type <- retained_registry$data_type[match(covariate_name, retained_registry$covariate)]
      association <- pairwise_covariate_association(
        x = pc_table[[pc_name]],
        y = pc_table[[covariate_name]],
        x_name = pc_name,
        y_name = covariate_name,
        x_type = "numeric",
        y_type = covariate_type
      )

      association$pc <- pc_name
      association$covariate <- covariate_name
      association$covariate_role <- retained_registry$role[match(covariate_name, retained_registry$covariate)]
      association$covariate_type <- covariate_type
      results[[result_index]] <- association
      result_index <- result_index + 1L
    }
  }

  do.call(rbind, results)
}

build_heatmap_matrix <- function(long_table, row_column, column_column, value_column, fill_value = NA_real_) {
  row_values <- unique(long_table[[row_column]])
  column_values <- unique(long_table[[column_column]])
  matrix_out <- matrix(fill_value, nrow = length(row_values), ncol = length(column_values))
  rownames(matrix_out) <- row_values
  colnames(matrix_out) <- column_values

  for (i in seq_len(nrow(long_table))) {
    row_name <- long_table[[row_column]][[i]]
    col_name <- long_table[[column_column]][[i]]
    matrix_out[row_name, col_name] <- long_table[[value_column]][[i]]
  }

  matrix_out
}

save_numeric_correlation_heatmap <- function(coldata, retained_registry, output_file) {
  numeric_covariates <- retained_registry$covariate[
    retained_registry$retain_for_pc_screening &
      retained_registry$data_type == "numeric"
  ]

  if (length(numeric_covariates) < 2) {
    return(invisible(NULL))
  }

  correlation_matrix <- stats::cor(
    as.data.frame(coldata[, numeric_covariates, drop = FALSE]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )

  grDevices::pdf(output_file, width = 9, height = 8)
  pheatmap::pheatmap(
    correlation_matrix,
    main = "Spearman correlations among retained numeric covariates",
    color = grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
    breaks = seq(-1, 1, length.out = 101)
  )
  grDevices::dev.off()
}

save_matrix_heatmap <- function(matrix_data, output_file, title_text, palette_values, breaks = NA) {
  plot_matrix <- matrix_data
  plot_matrix[is.na(plot_matrix)] <- 0

  grDevices::pdf(output_file, width = 10, height = 8)
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

save_pca_scree_plot <- function(percent_var, output_file) {
  scree_table <- data.frame(
    pc = paste0("PC", seq_along(percent_var)),
    variance_percent = percent_var,
    cumulative_variance_percent = cumsum(percent_var),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(scree_table, ggplot2::aes(x = pc, y = variance_percent, group = 1)) +
    ggplot2::geom_col(fill = "#3b6fb6") +
    ggplot2::geom_line(ggplot2::aes(y = cumulative_variance_percent), color = "#b22222") +
    ggplot2::geom_point(ggplot2::aes(y = cumulative_variance_percent), color = "#b22222") +
    exp383_theme(base_size = 12) +
    ggplot2::labs(
      title = "Blind PCA scree plot",
      x = "Principal component",
      y = "Variance explained (%)"
    )

  ggplot2::ggsave(output_file, plot = p, width = 8, height = 5)
}

run_univariate_variance_partition <- function(vsd, retained_registry, coldata, output_root) {
  if (!requireNamespace("variancePartition", quietly = TRUE)) {
    writeLines(
      "Skipped because package 'variancePartition' is not installed in the active R environment.",
      con = file.path(output_root, "variance_partition_status.txt")
    )
    return(invisible(NULL))
  }

  assay_matrix <- SummarizedExperiment::assay(vsd)

  # VariancePartition is informative here, but on the full filtered matrix it
  # becomes slow enough to make a design-scout script clumsy to iterate with.
  # For scouting we only need a stable ranking of major sources of variation, so
  # we restrict the fit to the most variable genes.
  max_genes_for_variance_partition <- 1000L
  genes_before_subset <- nrow(assay_matrix)
  if (nrow(assay_matrix) > max_genes_for_variance_partition) {
    row_variances <- apply(assay_matrix, 1, stats::var, na.rm = TRUE)
    keep_idx <- order(row_variances, decreasing = TRUE)[seq_len(max_genes_for_variance_partition)]
    assay_matrix <- assay_matrix[keep_idx, , drop = FALSE]
  }

  variance_partition_input_summary <- data.frame(
    genes_before_subset = genes_before_subset,
    genes_used_for_variance_partition = nrow(assay_matrix),
    stringsAsFactors = FALSE
  )
  write_tsv(
    variance_partition_input_summary,
    file.path(output_root, "variance_partition_input_summary.tsv")
  )

  # For the variancePartition scout, keep the primary biological factor plus the
  # retained technical covariates. We intentionally exclude the decomposed
  # biological terms (`inoculum`, `dpi_factor`) here because they are nested
  # inside `group_assignment` and just duplicate the same signal in a slower,
  # less interpretable way.
  retained_covariates <- retained_registry$covariate[
    retained_registry$retain_for_pc_screening &
      (retained_registry$role != "biological" | retained_registry$covariate == "group_assignment")
  ]

  results <- list()
  result_index <- 1L
  for (covariate_name in retained_covariates) {
    covariate_values <- coldata[[covariate_name]]
    non_missing_values <- stats::na.omit(covariate_values)

    if (length(unique(non_missing_values)) < 2) {
      next
    }

    vp_data <- data.frame(value = covariate_values)
    rownames(vp_data) <- rownames(coldata)

    vp_fit <- tryCatch(
      variancePartition::fitExtractVarPartModel(
        exprObj = assay_matrix,
        formula = stats::as.formula("~ value"),
        data = vp_data
      ),
      error = function(e) e
    )

    if (inherits(vp_fit, "error")) {
      results[[result_index]] <- data.frame(
        covariate = covariate_name,
        mean_fraction_variance = NA_real_,
        median_fraction_variance = NA_real_,
        p90_fraction_variance = NA_real_,
        status = paste("error:", conditionMessage(vp_fit)),
        stringsAsFactors = FALSE
      )
      result_index <- result_index + 1L
      next
    }

    value_column <- setdiff(colnames(vp_fit), "Residuals")
    value_column <- value_column[[1]]
    value_scores <- vp_fit[[value_column]]

    results[[result_index]] <- data.frame(
      covariate = covariate_name,
      mean_fraction_variance = mean(value_scores, na.rm = TRUE),
      median_fraction_variance = stats::median(value_scores, na.rm = TRUE),
      p90_fraction_variance = as.numeric(stats::quantile(value_scores, probs = 0.9, na.rm = TRUE)),
      status = "ok",
      stringsAsFactors = FALSE
    )
    result_index <- result_index + 1L
  }

  if (length(results) == 0) {
    writeLines(
      "Skipped because no retained covariates had sufficient variation for univariate variancePartition fits.",
      con = file.path(output_root, "variance_partition_status.txt")
    )
    return(invisible(NULL))
  }

  variance_partition_table <- do.call(rbind, results)
  variance_partition_table <- variance_partition_table[
    order(variance_partition_table$mean_fraction_variance, decreasing = TRUE, na.last = TRUE),
    ,
    drop = FALSE
  ]

  write_tsv(variance_partition_table, file.path(output_root, "variance_partition_univariate.tsv"))

  plot_table <- subset(variance_partition_table, status == "ok" & !is.na(mean_fraction_variance))
  if (nrow(plot_table) > 0) {
    p <- ggplot2::ggplot(
      plot_table,
      ggplot2::aes(
        x = reorder(covariate, mean_fraction_variance),
        y = mean_fraction_variance
      )
    ) +
      ggplot2::geom_col(fill = "#3b6fb6") +
      ggplot2::coord_flip() +
      exp383_theme(base_size = 11) +
      ggplot2::labs(
        title = "Univariate variancePartition summary",
        x = "Covariate",
        y = "Mean fraction variance explained"
      )

    ggplot2::ggsave(
      file.path(output_root, "variance_partition_univariate.pdf"),
      plot = p,
      width = 8,
      height = 6
    )
  }
}

# -----------------
# Main script logic
# -----------------
args <- parse_cli_args()
input_rds <- normalizePath(args$input_rds, winslash = "/", mustWork = FALSE)
stop_if_missing(input_rds, "Prepared input RDS")

design_table <- load_design_table(args)
analysis_input <- readRDS(input_rds)

if (!all(c("txi_gene", "sample_metadata", "gene_annotation") %in% names(analysis_input))) {
  stop(
    "Input RDS does not look like the output of 01_build_salmon_gene_inputs.R",
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

if (identical(args$min_samples, "auto")) {
  group_sizes <- table(coldata_pop[[args$filter_group_var]])
  min_samples <- as.integer(min(group_sizes))
} else {
  min_samples <- args$min_samples
}

keep_genes <- compute_shared_filter(
  count_matrix = txi_pop$counts,
  coldata = coldata_pop,
  min_count = args$min_count,
  min_samples = min_samples,
  filter_group_var = args$filter_group_var
)

txi_pop_filtered <- subset_txi(
  txi = txi_pop,
  keep_rows = keep_genes,
  keep_cols = rep(TRUE, ncol(txi_pop$counts))
)

gene_annotation <- analysis_input$gene_annotation
shared_root <- file.path(args$output_root, args$population, "shared")
dir.create(shared_root, recursive = TRUE, showWarnings = FALSE)

filter_summary <- data.frame(
  population = args$population,
  n_samples = nrow(coldata_pop),
  n_genes_before_filter = nrow(txi_pop$counts),
  n_genes_after_filter = nrow(txi_pop_filtered$counts),
  min_count = args$min_count,
  min_samples = min_samples,
  filter_group_var = args$filter_group_var,
  stringsAsFactors = FALSE
)

write_tsv(coldata_pop, file.path(shared_root, "population_metadata.tsv"))
write_tsv(filter_summary, file.path(shared_root, "filter_summary.tsv"))
write_tsv(
  data.frame(gene_id = rownames(txi_pop_filtered$counts), stringsAsFactors = FALSE),
  file.path(shared_root, "filtered_gene_ids.tsv")
)

# ----------------------------
# Shared covariate-screen step
# ----------------------------
message("Running shared covariate screening for population: ", args$population)

screen_root <- file.path(shared_root, "covariate_screening")
screen_plot_dir <- file.path(screen_root, "plots")
screen_table_dir <- file.path(screen_root, "tables")
dir.create(screen_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(screen_table_dir, recursive = TRUE, showWarnings = FALSE)

registry <- build_covariate_registry(coldata_pop)
pairwise_associations <- build_pairwise_covariate_associations(coldata_pop, registry)
numeric_collinearity_summary <- build_numeric_collinearity_summary(pairwise_associations, registry)
retained_covariates <- build_retained_covariate_table(registry, numeric_collinearity_summary)

dds_blind <- DESeq2::DESeqDataSetFromTximport(
  txi = txi_pop_filtered,
  colData = coldata_pop,
  design = stats::as.formula("~ 1")
)
vsd_blind <- DESeq2::vst(dds_blind, blind = TRUE)
blind_pca_table <- build_pc_score_table(vsd_blind, coldata_pop, n_pcs = 10)

pc_variance_table <- data.frame(
  pc = paste0("PC", seq_along(attr(blind_pca_table, "percent_var"))),
  variance_percent = round(attr(blind_pca_table, "percent_var"), 4),
  cumulative_variance_percent = round(cumsum(attr(blind_pca_table, "percent_var")), 4),
  stringsAsFactors = FALSE
)

pc_covariate_associations <- build_pc_covariate_associations(blind_pca_table, retained_covariates)

write_tsv(registry, file.path(screen_table_dir, "candidate_covariates.tsv"))
write_tsv(pairwise_associations, file.path(screen_table_dir, "pairwise_covariate_associations.tsv"))
write_tsv(numeric_collinearity_summary, file.path(screen_table_dir, "technical_numeric_collinearity.tsv"))
write_tsv(retained_covariates, file.path(screen_table_dir, "retained_covariates.tsv"))
write_tsv(blind_pca_table, file.path(screen_table_dir, "blind_pca_scores.tsv"))
write_tsv(pc_variance_table, file.path(screen_table_dir, "pc_variance.tsv"))
write_tsv(pc_covariate_associations, file.path(screen_table_dir, "pc_covariate_associations.tsv"))

save_pca_scree_plot(
  percent_var = attr(blind_pca_table, "percent_var"),
  output_file = file.path(screen_plot_dir, "blind_pca_scree.pdf")
)

save_numeric_correlation_heatmap(
  coldata = coldata_pop,
  retained_registry = retained_covariates,
  output_file = file.path(screen_plot_dir, "retained_numeric_covariate_correlations.pdf")
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
  output_file = file.path(screen_plot_dir, "pc_covariate_effect_sizes.pdf"),
  title_text = "PC-covariate association effect sizes",
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
  output_file = file.path(screen_plot_dir, "pc_covariate_neglog10_pvalues.pdf"),
  title_text = "PC-covariate association significance (-log10 p)",
  palette_values = grDevices::colorRampPalette(c("white", "#2166ac"))(100)
)

blind_pca_plot_table <- blind_pca_table
attr(blind_pca_plot_table, "percent_var") <- attr(blind_pca_table, "percent_var")[1:2]
shared_plot_covariates <- retained_covariates$covariate[retained_covariates$retain_for_pc_screening]
for (var_name in shared_plot_covariates) {
  save_pca_plot(
    pca_table = blind_pca_plot_table,
    color_var = var_name,
    title_text = sprintf("Blind PCA: %s", var_name),
    output_file = file.path(screen_plot_dir, sprintf("blind_pca_by_%s.pdf", var_name))
  )
}

run_univariate_variance_partition(
  vsd = vsd_blind,
  retained_registry = retained_covariates,
  coldata = coldata_pop,
  output_root = screen_root
)

# ---------------------------------------
# Design-specific DESeq2 QC and behaviour
# ---------------------------------------
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
      paste0("min_samples: ", min_samples),
      paste0("filter_group_var: ", args$filter_group_var)
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

  pca_table <- build_pca_table(vsd, coldata_pop)
  attr(pca_table, "percent_var") <- round(
    100 * attr(DESeq2::plotPCA(vsd, intgroup = "group_assignment", returnData = TRUE), "percentVar"),
    2
  )

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
      c("group_assignment", "inoculum", "dpi_factor", "inoculation_batch"),
      colnames(coldata_pop)
    ), drop = FALSE]),
    output_file = file.path(plot_dir, "sample_distance_heatmap.pdf")
  )
  save_dispersion_plot(dds, file.path(plot_dir, "dispersion_plot.pdf"))
  save_size_factor_plot(size_factor_table, file.path(plot_dir, "size_factors.pdf"))

  variables_to_plot <- unique(c(
    "group_assignment",
    all.vars(design_formula_obj),
    "inoculum",
    "dpi_factor",
    "inoculation_batch",
    "sample_mass_mg",
    "incubation_time_hrs",
    "date_nuc_prep_days",
    "rrna_alignment_percent",
    "salmon_percent_mapped",
    "fastp_percent_surviving"
  ))

  for (var_name in variables_to_plot) {
    if (!var_name %in% colnames(pca_table)) {
      next
    }

    save_pca_plot(
      pca_table = pca_table,
      color_var = var_name,
      title_text = sprintf("PCA: %s | %s", design_id, var_name),
      output_file = file.path(plot_dir, sprintf("pca_by_%s.pdf", var_name))
    )
  }

  writeLines(capture.output(sessionInfo()), con = file.path(design_dir, "sessionInfo.txt"))
}

message("Finished DESeq2 design scouting for population: ", args$population)
