# Utility functions for scripts/02b_covariate_screening_main_script.R.
# Keep executable workflow code in the main script; keep support code here.

# -----------------------------
# Command-line argument parsing
# -----------------------------
# Only the sorted population is exposed at the command line.
parse_cli_args <- function() {
  defaults <- list(
    population = NA_character_
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

prefix_plot_title <- function(population, title_text) {
  if (is.null(population) || is.na(population) || population == "") {
    return(title_text)
  }

  paste0(population, ": ", title_text)
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
set3_palette <- function(n) {
  base_palette <- c(
    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
    "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
  )

  if (n <= length(base_palette)) {
    return(base_palette[seq_len(n)])
  }

  grDevices::colorRampPalette(base_palette)(n)
}

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

  color_values <- pca_table[[color_var]]
  is_numeric_color <- is.numeric(color_values) || is.integer(color_values)
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

  if (is_numeric_color) {
    p <- p + ggplot2::scale_colour_viridis_c(option = "viridis", na.value = "grey70")
  } else if (identical(color_var, "group_assignment")) {
    color_levels <- order_known_levels(color_values, exp383_group_assignment_levels)
    p <- p + ggplot2::scale_colour_manual(
      values = exp383_group_assignment_palette[color_levels],
      breaks = color_levels,
      labels = format_group_assignment_label(color_levels),
      na.value = "grey70"
    )
  } else {
    color_levels <- if (is.factor(color_values)) {
      levels(droplevels(color_values))
    } else {
      unique(as.character(color_values[!is.na(color_values)]))
    }
    p <- p + ggplot2::scale_colour_manual(
      values = stats::setNames(set3_palette(length(color_levels)), color_levels),
      na.value = "grey70"
    )
  }

  exp383_save_ggplot(output_file, plot = p, width = 7, height = 5.5)
}

save_pca_scree_plot <- function(pc_variance_table, n_selected_pcs, variance_threshold, output_file, population = NULL) {
  threshold_percent <- 100 * variance_threshold
  p <- ggplot2::ggplot(pc_variance_table, ggplot2::aes(x = pc_index, y = variance_percent)) +
    ggplot2::geom_col(fill = "#3b6fb6") +
    ggplot2::geom_line(ggplot2::aes(y = cumulative_variance_percent), color = "#b22222", group = 1) +
    ggplot2::geom_point(ggplot2::aes(y = cumulative_variance_percent), color = "#b22222") +
    ggplot2::geom_vline(xintercept = n_selected_pcs + 0.5, linetype = 2, color = "grey40") +
    exp383_theme(base_size = 12) +
    ggplot2::scale_x_continuous(
      breaks = pc_variance_table$pc_index,
      guide = ggplot2::guide_axis(angle = 90)
    ) +
    ggplot2::labs(
      title = prefix_plot_title(population, "Blind PCA scree plot"),
      subtitle = sprintf(
        "Selected PCs explain at least %.0f%% cumulative variance (%d PCs retained)",
        threshold_percent,
        n_selected_pcs
      ),
      x = "Principal component",
      y = "Variance explained (%)"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1),
      plot.margin = ggplot2::margin(t = 8, r = 12, b = 12, l = 8)
    )

  exp383_save_ggplot(output_file, plot = p, width = 10, height = 5)
}

save_numeric_correlation_heatmap <- function(coldata, registry, output_file, population = NULL) {
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
    main = paste(
      prefix_plot_title(population, "Spearman correlations among retained technical numeric covariates"),
      "Variables shown after collinearity pruning.",
      sep = "\n"
    ),
    color = grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
    breaks = seq(-1, 1, length.out = 101)
  )
  grDevices::dev.off()
}

save_metric_barplot <- function(plot_table, metric_column, title_text, y_label, output_file, subtitle_text = NULL, population = NULL) {
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
      title = prefix_plot_title(population, title_text),
      subtitle = subtitle_text,
      x = "Covariate",
      y = y_label,
      fill = "Selected"
    )

  exp383_save_ggplot(output_file, plot = p, width = 8, height = 6)
}

save_variance_partition_plot <- function(varpart_object, output_file, population = NULL) {
  exp383_open_png_device(output_file, width = 10, height = 8)
  plot_object <- variancePartition::sortCols(varpart_object) |>
    variancePartition::plotVarPart(label.angle = 90)
  plot_object <- plot_object +
    ggplot2::labs(
      title = prefix_plot_title(population, "Multivariate variancePartition by gene"),
      subtitle = "Gene-level variance explained by each model term."
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.margin = ggplot2::margin(t = 12, r = 12, b = 32, l = 12)
    )
  print(plot_object)
  grDevices::dev.off()
}

# -------------------------------------------------
# Selected covariate balance against biology plots
# -------------------------------------------------
# These helpers check whether selected technical covariates are distributed
# unevenly across the biological variables used to define the DEA contrasts.
exp383_inoculum_levels <- c("RML", "ME7", "22L", "CBH")
exp383_dpi_levels <- c("60", "90", "120")

exp383_group_assignment_levels <- as.vector(outer(
  exp383_inoculum_levels,
  exp383_dpi_levels,
  FUN = function(inoculum, dpi) paste0("group_", inoculum, "_", dpi)
))

exp383_group_assignment_palette <- stats::setNames(
  set3_palette(length(exp383_group_assignment_levels)),
  exp383_group_assignment_levels
)

exp383_inoculum_palette <- c(
  "RML" = "#8DD3C7",
  "ME7" = "#FB8072",
  "22L" = "#B3DE69",
  "CBH" = "#FDB462"
)

exp383_dpi_palette <- c(
  "60" = "#6A3D9A",
  "90" = "#E31A1C",
  "120" = "#FF7F00"
)

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  gsub("_+", "_", x)
}

order_known_levels <- function(values, preferred_levels) {
  observed <- unique(as.character(values[!is.na(values)]))
  c(intersect(preferred_levels, observed), sort(setdiff(observed, preferred_levels)))
}

format_group_assignment_label <- function(x) {
  x <- sub("^group_", "", x)
  gsub("_", " ", x)
}

extract_group_assignment_parts <- function(group_assignment) {
  stripped <- sub("^group_", "", as.character(group_assignment))
  parts <- strsplit(stripped, "_", fixed = TRUE)
  data.frame(
    inoculum = vapply(parts, `[`, character(1), 1),
    dpi = vapply(parts, `[`, character(1), 2),
    stringsAsFactors = FALSE
  )
}

get_selected_technical_covariates <- function(selection_metrics) {
  selected <- selection_metrics[
    selection_metrics$is_candidate & selection_metrics$selected,
    c("covariate", "analysis_column", "data_type"),
    drop = FALSE
  ]
  selected[order(selected$data_type, selected$covariate), , drop = FALSE]
}

biological_covariate_table <- function(coldata) {
  data.frame(
    biological_variable = c("group_assignment", "dpi", "inoculum"),
    analysis_column = c("group_assignment", "dpi", "inoculum"),
    data_type = "categorical",
    stringsAsFactors = FALSE
  )
}

ordered_biological_values <- function(x, biological_variable) {
  if (biological_variable == "group_assignment") {
    levels <- order_known_levels(x, exp383_group_assignment_levels)
    return(factor(as.character(x), levels = levels))
  }

  if (biological_variable == "dpi") {
    levels <- order_known_levels(as.character(x), exp383_dpi_levels)
    return(factor(as.character(x), levels = levels))
  }

  if (biological_variable == "inoculum") {
    levels <- order_known_levels(x, exp383_inoculum_levels)
    return(factor(as.character(x), levels = levels))
  }

  factor(as.character(x))
}

build_selected_covariate_biology_associations <- function(coldata, selection_metrics) {
  selected_covariates <- get_selected_technical_covariates(selection_metrics)
  biological_covariates <- biological_covariate_table(coldata)

  empty_result <- data.frame(
    covariate = character(),
    analysis_column = character(),
    covariate_type = character(),
    biological_variable = character(),
    biological_column = character(),
    n_complete = integer(),
    test_name = character(),
    effect_size_name = character(),
    effect_size = numeric(),
    effect_size_abs = numeric(),
    statistic = numeric(),
    p_value = numeric(),
    padj = numeric(),
    stringsAsFactors = FALSE
  )

  if (nrow(selected_covariates) == 0) {
    return(empty_result)
  }

  results <- list()
  result_index <- 1L
  for (i in seq_len(nrow(selected_covariates))) {
    covariate_row <- selected_covariates[i, , drop = FALSE]

    for (j in seq_len(nrow(biological_covariates))) {
      biological_row <- biological_covariates[j, , drop = FALSE]
      biological_values <- ordered_biological_values(
        coldata[[biological_row$analysis_column]],
        biological_row$biological_variable
      )

      association <- pairwise_covariate_association(
        x = coldata[[covariate_row$analysis_column]],
        y = biological_values,
        x_name = covariate_row$covariate,
        y_name = biological_row$biological_variable,
        x_type = covariate_row$data_type,
        y_type = biological_row$data_type
      )

      results[[result_index]] <- data.frame(
        covariate = covariate_row$covariate,
        analysis_column = covariate_row$analysis_column,
        covariate_type = covariate_row$data_type,
        biological_variable = biological_row$biological_variable,
        biological_column = biological_row$analysis_column,
        n_complete = association$n_complete,
        test_name = association$test_name,
        effect_size_name = association$effect_size_name,
        effect_size = association$effect_size,
        effect_size_abs = association$effect_size_abs,
        statistic = association$statistic,
        p_value = association$p_value,
        padj = NA_real_,
        stringsAsFactors = FALSE
      )
      result_index <- result_index + 1L
    }
  }

  association_table <- do.call(rbind, results)
  association_table$padj <- stats::p.adjust(association_table$p_value, method = "fdr")
  association_table
}

build_selected_numeric_covariate_biology_data <- function(coldata, selection_metrics) {
  selected_covariates <- get_selected_technical_covariates(selection_metrics)
  selected_covariates <- selected_covariates[selected_covariates$data_type == "numeric", , drop = FALSE]

  empty_result <- data.frame(
    sample = character(),
    covariate = character(),
    analysis_column = character(),
    covariate_value = numeric(),
    group_assignment = character(),
    dpi = character(),
    inoculum = character(),
    stringsAsFactors = FALSE
  )

  if (nrow(selected_covariates) == 0) {
    return(empty_result)
  }

  sample_ids <- if ("sample" %in% colnames(coldata)) coldata$sample else rownames(coldata)
  results <- lapply(seq_len(nrow(selected_covariates)), function(i) {
    covariate_row <- selected_covariates[i, , drop = FALSE]
    data.frame(
      sample = as.character(sample_ids),
      covariate = covariate_row$covariate,
      analysis_column = covariate_row$analysis_column,
      covariate_value = as.numeric(coldata[[covariate_row$analysis_column]]),
      group_assignment = as.character(coldata$group_assignment),
      dpi = as.character(coldata$dpi),
      inoculum = as.character(coldata$inoculum),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

build_selected_categorical_covariate_biology_data <- function(coldata, selection_metrics) {
  selected_covariates <- get_selected_technical_covariates(selection_metrics)
  selected_covariates <- selected_covariates[selected_covariates$data_type == "categorical", , drop = FALSE]
  biological_covariates <- biological_covariate_table(coldata)

  empty_result <- data.frame(
    covariate = character(),
    analysis_column = character(),
    biological_variable = character(),
    biological_value = character(),
    covariate_value = character(),
    n_samples = integer(),
    biological_group_n = integer(),
    proportion_within_biological_group = numeric(),
    stringsAsFactors = FALSE
  )

  if (nrow(selected_covariates) == 0) {
    return(empty_result)
  }

  results <- list()
  result_index <- 1L
  for (i in seq_len(nrow(selected_covariates))) {
    covariate_row <- selected_covariates[i, , drop = FALSE]

    for (j in seq_len(nrow(biological_covariates))) {
      biological_row <- biological_covariates[j, , drop = FALSE]
      biological_values <- ordered_biological_values(
        coldata[[biological_row$analysis_column]],
        biological_row$biological_variable
      )
      covariate_values <- factor(as.character(coldata[[covariate_row$analysis_column]]))

      count_table <- as.data.frame.matrix(table(biological_values, covariate_values))
      biological_group_n <- rowSums(count_table)

      for (biological_value in rownames(count_table)) {
        for (covariate_value in colnames(count_table)) {
          n_samples <- count_table[biological_value, covariate_value]
          group_n <- biological_group_n[[biological_value]]
          results[[result_index]] <- data.frame(
            covariate = covariate_row$covariate,
            analysis_column = covariate_row$analysis_column,
            biological_variable = biological_row$biological_variable,
            biological_value = biological_value,
            covariate_value = covariate_value,
            n_samples = as.integer(n_samples),
            biological_group_n = as.integer(group_n),
            proportion_within_biological_group = if (group_n > 0) n_samples / group_n else NA_real_,
            stringsAsFactors = FALSE
          )
          result_index <- result_index + 1L
        }
      }
    }
  }

  do.call(rbind, results)
}

save_selected_numeric_covariate_biology_plot <- function(plot_data, biological_variable, output_file, population = NULL) {
  if (nrow(plot_data) == 0) {
    return(invisible(NULL))
  }

  plot_data <- plot_data[!is.na(plot_data$covariate_value), , drop = FALSE]
  if (nrow(plot_data) == 0) {
    return(invisible(NULL))
  }

  plot_data$biological_value <- ordered_biological_values(
    plot_data[[biological_variable]],
    biological_variable
  )

  x_label <- switch(
    biological_variable,
    group_assignment = "Group assignment",
    dpi = "DPI",
    inoculum = "Inoculum",
    biological_variable
  )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = biological_value, y = covariate_value, fill = biological_value)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, colour = "grey35", alpha = 0.82, show.legend = FALSE) +
    ggplot2::geom_jitter(width = 0.16, height = 0, size = 1.6, alpha = 0.72, colour = "#1f1f1f") +
    ggplot2::facet_wrap(~covariate, scales = "free_y", ncol = 2) +
    exp383_theme(base_size = 11) +
    ggplot2::labs(
      title = prefix_plot_title(population, sprintf("Selected numeric technical covariates by %s", x_label)),
      subtitle = sprintf("Z-scored selected technical covariates grouped by %s.", x_label),
      x = x_label,
      y = "Z-scored covariate value"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.margin = ggplot2::margin(t = 8, r = 12, b = 18, l = 8)
    )

  fill_levels <- levels(plot_data$biological_value)
  fill_values <- switch(
    biological_variable,
    group_assignment = exp383_group_assignment_palette[fill_levels],
    dpi = exp383_dpi_palette[fill_levels],
    inoculum = exp383_inoculum_palette[fill_levels],
    stats::setNames(set3_palette(length(fill_levels)), fill_levels)
  )

  p <- p + ggplot2::scale_fill_manual(values = fill_values, drop = FALSE, guide = "none")

  plot_width <- if (biological_variable == "group_assignment") 11 else 8
  if (biological_variable == "group_assignment") {
    p <- p + ggplot2::scale_x_discrete(labels = format_group_assignment_label)
  }

  exp383_save_ggplot(output_file, plot = p, width = plot_width, height = 6.5)
}

save_selected_categorical_covariate_group_heatmaps <- function(count_data, output_dir, population = NULL) {
  group_data <- count_data[count_data$biological_variable == "group_assignment", , drop = FALSE]
  if (nrow(group_data) == 0) {
    return(invisible(NULL))
  }

  covariates <- unique(group_data$covariate)
  for (covariate_name in covariates) {
    covariate_data <- group_data[group_data$covariate == covariate_name, , drop = FALSE]
    row_order <- order_known_levels(covariate_data$biological_value, exp383_group_assignment_levels)
    column_order <- sort(unique(covariate_data$covariate_value))

    heatmap_matrix <- matrix(
      NA_real_,
      nrow = length(row_order),
      ncol = length(column_order),
      dimnames = list(row_order, column_order)
    )

    for (i in seq_len(nrow(covariate_data))) {
      heatmap_matrix[
        covariate_data$biological_value[[i]],
        covariate_data$covariate_value[[i]]
      ] <- covariate_data$proportion_within_biological_group[[i]]
    }

    row_parts <- extract_group_assignment_parts(rownames(heatmap_matrix))
    row_annotation <- data.frame(
      inoculum = factor(row_parts$inoculum, levels = exp383_inoculum_levels),
      dpi = factor(row_parts$dpi, levels = exp383_dpi_levels),
      stringsAsFactors = FALSE
    )
    rownames(row_annotation) <- rownames(heatmap_matrix)

    # Match the proteomics EDA heatmap convention: fixed row order, no
    # clustering, visible group gaps, fixed annotation colours, and black NA.
    dpi_group_sizes <- table(row_annotation$dpi)
    row_gaps <- cumsum(as.integer(dpi_group_sizes))
    row_gaps <- row_gaps[row_gaps < nrow(heatmap_matrix)]

    annotation_colors <- list(
      inoculum = exp383_inoculum_palette[exp383_inoculum_levels],
      dpi = exp383_dpi_palette[exp383_dpi_levels]
    )

    display_matrix <- heatmap_matrix
    rownames(display_matrix) <- format_group_assignment_label(rownames(display_matrix))
    rownames(row_annotation) <- rownames(display_matrix)

    pheat <- pheatmap::pheatmap(
      display_matrix,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_colnames = TRUE,
      show_rownames = TRUE,
      fontsize_row = 8,
      fontsize_col = 9,
      main = paste(
        prefix_plot_title(population, sprintf("%s by group assignment", covariate_name)),
        "Cell values are row proportions within group assignment.",
        sep = "\n"
      ),
      annotation_row = row_annotation,
      annotation_colors = annotation_colors,
      annotation_names_row = FALSE,
      gaps_row = row_gaps,
      na_col = "black",
      silent = TRUE
    )

    output_file <- file.path(
      output_dir,
      sprintf("selected_categorical_covariate_by_group_assignment_%s.png", sanitize_filename(covariate_name))
    )
    exp383_open_png_device(output_file, width = 9, height = 7)
    grid::grid.newpage()
    grid::grid.draw(pheat$gtable)
    grDevices::dev.off()
  }
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
