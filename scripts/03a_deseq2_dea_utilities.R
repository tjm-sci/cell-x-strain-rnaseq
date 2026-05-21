# Utility functions for scripts/03b_deseq2_dea_fit_models.R and
# scripts/03c_deseq2_dea_export_results.R.
#
# Keep the production scripts short: path handling, filtering, DESeq2 model QC,
# contrast table handling, and plotting helpers live here.

# -----------------
# Small utilities
# -----------------
stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
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

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  gsub("_+", "_", x)
}

parse_key_value_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %% 2 != 0) {
    stop("Arguments must be provided as --key value pairs.", call. = FALSE)
  }

  parsed <- defaults
  if (length(args) == 0) {
    return(parsed)
  }

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

  parsed
}

parse_population_selection <- function(population_arg, available_populations) {
  available_populations <- as.character(available_populations)
  if (identical(population_arg, "all")) {
    return(available_populations)
  }

  requested <- strsplit(population_arg, ",", fixed = TRUE)[[1]]
  requested <- trimws(requested)
  missing <- setdiff(requested, available_populations)
  if (length(missing) > 0) {
    stop(
      sprintf("Requested population(s) not found: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  requested
}

# ----------------------------
# Shared biological aesthetics
# ----------------------------
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

format_group_assignment_label <- function(x) {
  x <- sub("^group_", "", as.character(x))
  gsub("_", " ", x)
}

order_known_levels <- function(values, preferred_levels) {
  observed <- unique(as.character(values[!is.na(values)]))
  c(intersect(preferred_levels, observed), sort(setdiff(observed, preferred_levels)))
}

# -----------------------------
# Metadata and tximport helpers
# -----------------------------
load_dea_input <- function(input_rds) {
  stop_if_missing(input_rds, "Corrected input RDS")
  analysis_input <- readRDS(input_rds)

  required_objects <- c("txi_gene", "sample_metadata", "gene_annotation")
  missing_objects <- setdiff(required_objects, names(analysis_input))
  if (length(missing_objects) > 0) {
    stop(
      "Input RDS is missing required object(s): ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  analysis_input
}

prepare_coldata_for_deseq <- function(coldata) {
  coldata <- as.data.frame(coldata, stringsAsFactors = FALSE)
  if ("sample" %in% colnames(coldata)) {
    rownames(coldata) <- coldata$sample
  }

  character_columns <- vapply(coldata, is.character, logical(1))
  coldata[character_columns] <- lapply(coldata[character_columns], factor)

  if ("population" %in% colnames(coldata)) {
    coldata$population <- factor(as.character(coldata$population), levels = c("NeuN", "SOX10", "SOX2", "PU1"))
  }
  if ("group_assignment" %in% colnames(coldata)) {
    group_levels <- order_known_levels(coldata$group_assignment, exp383_group_assignment_levels)
    coldata$group_assignment <- factor(as.character(coldata$group_assignment), levels = group_levels)
  }

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

ensure_formula_z_columns <- function(coldata, design_formula) {
  formula_terms <- all.vars(stats::as.formula(design_formula))
  z_terms <- formula_terms[grepl("_z$", formula_terms)]

  for (z_col in z_terms) {
    source_col <- sub("_z$", "", z_col)
    if (!z_col %in% colnames(coldata) && source_col %in% colnames(coldata)) {
      coldata[[z_col]] <- scale_numeric_covariate(coldata[[source_col]])
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

compute_group_count_filter <- function(count_matrix, grouping, min_count, min_samples) {
  split_indices <- split(seq_len(ncol(count_matrix)), grouping)

  keep <- rep(FALSE, nrow(count_matrix))
  for (group_name in names(split_indices)) {
    group_idx <- split_indices[[group_name]]
    keep_in_group <- rowSums(count_matrix[, group_idx, drop = FALSE] >= min_count) >= min_samples
    keep <- keep | keep_in_group
  }

  keep
}

read_population_design_formula <- function(screening_root, population) {
  formula_file <- file.path(
    screening_root,
    population,
    "covariate_screening",
    "recommended_design_formula.txt"
  )
  stop_if_missing(formula_file, sprintf("%s recommended design formula", population))

  formula_text <- trimws(readLines(formula_file, warn = FALSE))
  formula_text <- formula_text[nzchar(formula_text)]
  if (length(formula_text) != 1) {
    stop(sprintf("Expected exactly one design formula in %s", formula_file), call. = FALSE)
  }

  formula_text
}

prepare_population_for_deseq <- function(analysis_input, population, design_formula, min_count, min_samples) {
  coldata_all <- prepare_coldata_for_deseq(analysis_input$sample_metadata)
  coldata_all <- ensure_formula_z_columns(coldata_all, design_formula)

  if (!population %in% as.character(coldata_all$population)) {
    stop(sprintf("Population not found in metadata: %s", population), call. = FALSE)
  }

  keep_samples <- coldata_all$population == population
  coldata_pop <- droplevels(coldata_all[keep_samples, , drop = FALSE])
  rownames(coldata_pop) <- coldata_pop$sample

  required_terms <- all.vars(stats::as.formula(design_formula))
  missing_terms <- setdiff(required_terms, colnames(coldata_pop))
  if (length(missing_terms) > 0) {
    stop(
      sprintf("Design formula term(s) missing from %s metadata: %s", population, paste(missing_terms, collapse = ", ")),
      call. = FALSE
    )
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
    stringsAsFactors = FALSE
  )

  list(
    coldata = coldata_pop,
    txi = txi_filtered,
    keep_genes = keep_genes,
    filter_summary = filter_summary
  )
}

check_design_matrix_full_rank <- function(design_formula, coldata) {
  design_matrix <- stats::model.matrix(stats::as.formula(design_formula), data = as.data.frame(coldata))
  matrix_rank <- qr(design_matrix)$rank

  list(
    design_matrix = design_matrix,
    is_full_rank = matrix_rank == ncol(design_matrix),
    rank = matrix_rank,
    n_columns = ncol(design_matrix)
  )
}

get_sample_scaling_factors <- function(dds) {
  size_factors <- DESeq2::sizeFactors(dds)
  if (!is.null(size_factors)) {
    return(size_factors)
  }

  normalization_factors <- DESeq2::normalizationFactors(dds)
  if (is.null(normalization_factors)) {
    stop("DESeq2 object has neither sizeFactors nor normalizationFactors.", call. = FALSE)
  }

  exp(colMeans(log(normalization_factors)))
}

# ------------------------
# Model QC plot functions
# ------------------------
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

save_model_pca_plot <- function(pca_table, color_var, title_text, output_file) {
  if (!color_var %in% colnames(pca_table)) {
    return(invisible(NULL))
  }

  color_values <- pca_table[[color_var]]
  percent_var <- attr(pca_table, "percent_var")

  p <- ggplot2::ggplot(
    pca_table,
    ggplot2::aes(x = PC1, y = PC2, color = .data[[color_var]])
  ) +
    ggplot2::geom_point(size = 2.8, alpha = 0.9) +
    exp383_theme(base_size = 12) +
    ggplot2::labs(
      title = title_text,
      subtitle = "VST calculated from the fitted DESeq2 model.",
      x = sprintf("PC1 (%.2f%% variance)", percent_var[[1]]),
      y = sprintf("PC2 (%.2f%% variance)", percent_var[[2]]),
      color = color_var
    )

  if (is.numeric(color_values) || is.integer(color_values)) {
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

save_sample_distance_heatmap <- function(vsd, coldata, output_file, population) {
  assay_matrix <- SummarizedExperiment::assay(vsd)
  sample_dist <- as.matrix(stats::dist(t(assay_matrix)))

  sample_order <- order(coldata$group_assignment, coldata$sample)
  sample_ids <- rownames(coldata)[sample_order]
  sample_dist <- sample_dist[sample_ids, sample_ids, drop = FALSE]
  diag(sample_dist) <- NA

  annotation_df <- data.frame(
    group_assignment = coldata$group_assignment[sample_order],
    row.names = sample_ids
  )

  annotation_colors <- list(
    group_assignment = exp383_group_assignment_palette[levels(droplevels(annotation_df$group_assignment))]
  )

  group_gaps <- cumsum(as.integer(table(annotation_df$group_assignment)))
  group_gaps <- group_gaps[group_gaps < nrow(sample_dist)]

  pheat <- pheatmap::pheatmap(
    sample_dist,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 6,
    main = paste0(population, ": VST sample distances"),
    annotation_row = annotation_df,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    gaps_row = group_gaps,
    gaps_col = group_gaps,
    na_col = "black",
    silent = TRUE
  )

  exp383_open_png_device(output_file, width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(pheat$gtable)
  grDevices::dev.off()
}

save_size_factor_plot <- function(size_factor_table, output_file, population) {
  p <- ggplot2::ggplot(
    size_factor_table,
    ggplot2::aes(x = reorder(sample, size_factor), y = size_factor, fill = group_assignment)
  ) +
    ggplot2::geom_col(show.legend = TRUE) +
    ggplot2::coord_flip() +
    exp383_theme(base_size = 10) +
    ggplot2::scale_fill_manual(
      values = exp383_group_assignment_palette[levels(droplevels(size_factor_table$group_assignment))],
      labels = format_group_assignment_label
    ) +
    ggplot2::labs(
      title = paste0(population, ": DESeq2 size factors"),
      x = "Sample",
      y = "Size factor",
      fill = "Group"
    )

  exp383_save_ggplot(output_file, plot = p, width = 8, height = 10)
}

save_dispersion_plot <- function(dds, output_file, population) {
  exp383_open_png_device(output_file, width = 7, height = 5.5)
  DESeq2::plotDispEsts(dds, main = paste0(population, ": dispersion estimates"))
  grDevices::dev.off()
}

# ------------------------
# Contrast result helpers
# ------------------------
read_contrast_table <- function(contrast_tsv) {
  stop_if_missing(contrast_tsv, "DEA contrast TSV")
  contrast_table <- utils::read.delim(contrast_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  required_columns <- c("contrast_id", "contrast_family", "numerator_group", "denominator_group", "contrast_label")
  missing_columns <- setdiff(required_columns, colnames(contrast_table))
  if (length(missing_columns) > 0) {
    stop(
      sprintf("Contrast TSV missing required columns: %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  contrast_table$contrast_id <- vapply(contrast_table$contrast_id, sanitize_filename, character(1))
  contrast_table
}

annotate_deseq_results <- function(result_table, gene_annotation) {
  result_table <- as.data.frame(result_table, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"gene_id" %in% colnames(result_table)) {
    result_table$gene_id <- rownames(result_table)
  }

  annotation <- as.data.frame(gene_annotation, stringsAsFactors = FALSE, check.names = FALSE)
  annotation <- dplyr::transmute(
    annotation,
    gene_id = .data$gene_id,
    gene_symbol = .data$gene_name
  )

  result_table <- dplyr::left_join(result_table, annotation, by = "gene_id")
  result_table$gene_symbol <- ifelse(
    is.na(result_table$gene_symbol) | result_table$gene_symbol == "",
    result_table$gene_id,
    result_table$gene_symbol
  )

  result_table |>
    dplyr::select(
      "gene_symbol",
      "gene_id",
      dplyr::everything(),
      -dplyr::any_of(c("gene_name", "gene_label"))
    )
}

classify_de_direction <- function(log2_fold_change, padj, padj_cutoff) {
  ifelse(
    !is.na(padj) & padj <= padj_cutoff & log2_fold_change > 0,
    "up",
    ifelse(
      !is.na(padj) & padj <= padj_cutoff & log2_fold_change < 0,
      "down",
      "not_significant"
    )
  )
}

save_volcano_plot <- function(result_table, contrast_label, output_file, padj_cutoff, lfc_threshold, label_top_n,
                              p_mode = c("adjusted", "raw")) {
  p_mode <- match.arg(p_mode)
  plot_df <- result_table

  # The result table should expose pvalue/padj from the test used for calling
  # genes. For this analysis that is the DESeq2 non-zero LFC threshold test.
  plot_df$direction <- classify_de_direction(
    plot_df$log2FoldChange,
    plot_df$padj,
    padj_cutoff = padj_cutoff
  )

  if (p_mode == "adjusted") {
    plot_df$p_for_plot <- ifelse(
      is.na(plot_df$padj),
      NA_real_,
      pmax(plot_df$padj, .Machine$double.xmin)
    )
    y_axis_label <- "-log10(adjusted p-value)"
    p_line_values <- data.frame(
      threshold = -log10(padj_cutoff),
      line_type = "Threshold-test FDR = 0.05",
      stringsAsFactors = FALSE
    )
    subtitle_text <- sprintf("DESeq2 threshold test: |log2FC| > %.2f; FDR < %.2f", lfc_threshold, padj_cutoff)
  } else {
    plot_df$p_for_plot <- ifelse(
      is.na(plot_df$pvalue),
      NA_real_,
      pmax(plot_df$pvalue, .Machine$double.xmin)
    )

    fdr_hit_pvalues <- plot_df$pvalue[
      !is.na(plot_df$padj) &
        plot_df$padj <= padj_cutoff &
        !is.na(plot_df$pvalue)
    ]
    bh_raw_p_cutoff <- if (length(fdr_hit_pvalues) > 0) {
      max(fdr_hit_pvalues, na.rm = TRUE)
    } else {
      NA_real_
    }

    p_line_values <- data.frame(
      threshold = -log10(0.05),
      line_type = "Raw p = 0.05",
      stringsAsFactors = FALSE
    )
    if (!is.na(bh_raw_p_cutoff) && is.finite(bh_raw_p_cutoff)) {
      p_line_values <- rbind(
        p_line_values,
        data.frame(
          threshold = -log10(max(bh_raw_p_cutoff, .Machine$double.xmin)),
          line_type = "BH raw-p cutoff",
          stringsAsFactors = FALSE
        )
      )
    }

    y_axis_label <- "-log10(raw p-value)"
    subtitle_text <- "DESeq2 threshold-test raw p-values; green marks largest raw p among FDR-significant genes."
  }

  plot_df$neg_log10_p <- -log10(plot_df$p_for_plot)
  plot_df <- plot_df[
    !is.na(plot_df$log2FoldChange) &
      !is.na(plot_df$neg_log10_p),
    ,
    drop = FALSE
  ]

  label_df <- plot_df[
    plot_df$direction != "not_significant" &
      !is.na(plot_df$padj) &
      !is.na(plot_df$log2FoldChange),
    ,
    drop = FALSE
  ]
  label_df <- label_df[order(label_df$padj, -abs(label_df$log2FoldChange)), , drop = FALSE]
  label_df <- head(label_df, label_top_n)

  max_abs_lfc <- max(abs(plot_df$log2FoldChange), na.rm = TRUE)
  max_abs_lfc <- ifelse(is.finite(max_abs_lfc), max(max_abs_lfc, lfc_threshold), lfc_threshold)
  max_y <- max(c(plot_df$neg_log10_p, p_line_values$threshold), na.rm = TRUE)
  max_y <- ifelse(is.finite(max_y), max_y, 1)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = log2FoldChange, y = neg_log10_p)) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "grey40") +
    ggplot2::geom_hline(
      data = p_line_values,
      ggplot2::aes(yintercept = threshold, color = line_type),
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(ggplot2::aes(color = direction), size = 1.2, alpha = 0.8) +
    ggrepel::geom_label_repel(
      data = label_df,
      ggplot2::aes(label = gene_symbol),
      size = 3,
      box.padding = 0.35,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    ggplot2::coord_cartesian(xlim = c(-max_abs_lfc, max_abs_lfc), ylim = c(0, max_y)) +
    exp383_theme(base_size = 12) +
    ggplot2::scale_color_manual(
      values = c(
        "down" = "blue",
        "not_significant" = "grey70",
        "up" = "red",
        "Threshold-test FDR = 0.05" = "grey40",
        "Raw p = 0.05" = "black",
        "BH raw-p cutoff" = "#2CA25F"
      ),
      breaks = c("up", "down", "not_significant", "Threshold-test FDR = 0.05", "Raw p = 0.05", "BH raw-p cutoff")
    ) +
    ggplot2::labs(
      title = contrast_label,
      subtitle = subtitle_text,
      x = "log2 fold change",
      y = y_axis_label,
      color = NULL
    )

  exp383_save_ggplot(output_file, plot = p, width = 7, height = 5.5)
}

save_ma_plot <- function(result_table, contrast_label, output_file, padj_cutoff, lfc_threshold) {
  plot_df <- result_table
  plot_df$baseMean_plot <- ifelse(
    is.na(plot_df$baseMean),
    NA_real_,
    pmax(plot_df$baseMean, .Machine$double.xmin)
  )
  plot_df$direction <- classify_de_direction(
    plot_df$log2FoldChange,
    plot_df$padj,
    padj_cutoff = padj_cutoff
  )
  plot_df <- plot_df[
    !is.na(plot_df$baseMean_plot) &
      !is.na(plot_df$log2FoldChange),
    ,
    drop = FALSE
  ]

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = baseMean_plot, y = log2FoldChange)) +
    ggplot2::geom_hline(yintercept = c(-lfc_threshold, 0, lfc_threshold), linetype = "dashed", color = "grey55") +
    ggplot2::geom_point(ggplot2::aes(color = direction), size = 1.1, alpha = 0.75) +
    ggplot2::scale_x_log10() +
    exp383_theme(base_size = 12) +
    ggplot2::scale_color_manual(
      values = c("down" = "blue", "not_significant" = "grey70", "up" = "red"),
      breaks = c("up", "down", "not_significant")
    ) +
    ggplot2::labs(
      title = paste0(contrast_label, " MA plot"),
      x = "Mean normalised count",
      y = "log2 fold change",
      color = "Direction"
    )

  exp383_save_ggplot(output_file, plot = p, width = 7, height = 5.5)
}
