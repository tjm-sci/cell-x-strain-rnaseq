#!/usr/bin/env Rscript

# This script builds a whole-dataset expression overview for all 286 samples.
#
# It is intentionally separate from the population-specific DEA workflow. The
# goal here is not differential testing, but a project-level QC view that helps
# answer questions such as:
# - do the four FANS-derived nuclei populations separate cleanly?
# - are there obvious outliers or mis-sorted samples?
# - do key technical variables line up with the dominant axes of variation?
#
# The workflow is deliberately simple:
# 1. load the unified Salmon / tximport handoff object from script 01
# 2. apply one light global expression filter
# 3. build a DESeq2 object with design ~ 1
# 4. run a blind variance stabilizing transform
# 5. make PCA, sample-distance, and marker-gene QC plots
#
# No DEA is performed here.

suppressPackageStartupMessages({
  required_packages <- c(
    "DESeq2",
    "dplyr",
    "ggplot2",
    "pheatmap",
    "stringr",
    "tibble",
    "tidyr"
  )

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

# Load the shared project plotting helpers from the same scripts directory.
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
    output_dir = "results/dea/01b_global_sample_overview",
    min_count = "10",
    min_samples = "6",
    filter_group_var = "population",
    pca_ntop = "500"
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

  parsed$min_count <- as.numeric(parsed$min_count)
  parsed$min_samples <- as.integer(parsed$min_samples)
  parsed$pca_ntop <- as.integer(parsed$pca_ntop)

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

metadata_to_tibble <- function(coldata) {
  metadata_tbl <- tibble::as_tibble(coldata)

  # The handoff object already carries a `sample` column. We keep it as the
  # canonical sample identifier and only fall back to row names if an older
  # object lacks that column.
  if (!"sample" %in% colnames(metadata_tbl)) {
    metadata_tbl <- tibble::rownames_to_column(metadata_tbl, "sample")
  }

  metadata_tbl
}

subset_txi <- function(txi, keep_rows, keep_cols) {
  out <- txi
  out$counts <- txi$counts[keep_rows, keep_cols, drop = FALSE]
  out$abundance <- txi$abundance[keep_rows, keep_cols, drop = FALSE]
  out$length <- txi$length[keep_rows, keep_cols, drop = FALSE]
  out$countsFromAbundance <- txi$countsFromAbundance
  out
}

# The whole-dataset overview should still keep population-specific genes.
# Filtering within population levels avoids dropping genes that are strongly
# expressed in one sorted population but quiet in the others.
compute_grouped_filter <- function(count_matrix, grouping, min_count, min_samples) {
  split_indices <- split(seq_len(ncol(count_matrix)), grouping)

  keep <- rep(FALSE, nrow(count_matrix))
  for (group_name in names(split_indices)) {
    group_idx <- split_indices[[group_name]]
    keep_in_group <- rowSums(count_matrix[, group_idx, drop = FALSE] >= min_count) >= min_samples
    keep <- keep | keep_in_group
  }

  keep
}

# DESeq2's plotPCA() uses the most variable genes. We mimic that logic here so
# the overview behaves like a standard DESeq2 PCA, while still returning a tidy
# table that can be reused for multiple plots.
build_pca_coordinates <- function(vsd, coldata, ntop) {
  transformed_matrix <- SummarizedExperiment::assay(vsd)
  row_variances <- apply(transformed_matrix, 1, stats::var, na.rm = TRUE)
  keep_idx <- order(row_variances, decreasing = TRUE)[seq_len(min(ntop, nrow(transformed_matrix)))]

  pca <- stats::prcomp(t(transformed_matrix[keep_idx, , drop = FALSE]), center = TRUE, scale. = FALSE)
  percent_variance <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

  pca_table <- as.data.frame(pca$x[, seq_len(min(10, ncol(pca$x))), drop = FALSE]) |>
    tibble::rownames_to_column("sample") |>
    dplyr::left_join(
      metadata_to_tibble(coldata),
      by = "sample"
    )

  attr(pca_table, "percent_variance") <- percent_variance
  pca_table
}

save_pca_plot <- function(pca_table, color_var, percent_variance, output_file, ntop) {
  if (!color_var %in% colnames(pca_table)) {
    return(invisible(NULL))
  }

  is_numeric_colour <- is.numeric(pca_table[[color_var]])

  p <- ggplot2::ggplot(
    pca_table,
    ggplot2::aes(x = .data$PC1, y = .data$PC2, colour = .data[[color_var]])
  ) +
    ggplot2::geom_point(size = 2.3, alpha = 0.85) +
    exp383_theme(base_size = 12) +
    ggplot2::labs(
      title = paste("Global PCA coloured by", color_var),
      subtitle = sprintf("PCs calculated using the top %d most variable genes", ntop),
      x = sprintf("PC1 (%.1f%% variance)", percent_variance[[1]]),
      y = sprintf("PC2 (%.1f%% variance)", percent_variance[[2]]),
      colour = color_var
    )

  if (identical(color_var, "population")) {
    p <- p + exp383_scale_colour_population(name = "population")
  } else if (is_numeric_colour) {
    p <- p + ggplot2::scale_colour_gradient(low = "#d9e6f5", high = "#08306b")
  }

  ggplot2::ggsave(output_file, plot = p, width = 7.5, height = 6)
}

build_distance_legend_spec <- function(distance_matrix) {
  off_diagonal_values <- distance_matrix[upper.tri(distance_matrix)]
  if (length(off_diagonal_values) == 0) {
    off_diagonal_values <- 0
  }

  distance_max <- max(off_diagonal_values, na.rm = TRUE)

  list(
    breaks = c(0, distance_max),
    labels = c(
      "0",
      format(round(distance_max, 1), trim = TRUE, nsmall = 1)
    )
  )
}

save_sample_distance_heatmap <- function(vsd, coldata, output_file) {
  distance_matrix <- as.matrix(stats::dist(t(SummarizedExperiment::assay(vsd))))
  rownames(distance_matrix) <- rownames(coldata)
  colnames(distance_matrix) <- rownames(coldata)
  distance_legend <- build_distance_legend_spec(distance_matrix)

  # This script is only meant to show population-level sorting success. Keeping
  # just the population annotation avoids unreadable legends from many-level
  # experimental variables such as group assignment or inoculation batch.
  annotation_df <- coldata |>
    dplyr::transmute(population = exp383_population_display_factor(.data$population)) |>
    as.data.frame()

  grDevices::pdf(output_file, width = 11, height = 10)
  pheatmap::pheatmap(
    distance_matrix,
    annotation_col = annotation_df,
    annotation_row = annotation_df,
    annotation_colors = exp383_population_annotation_colors(),
    show_rownames = FALSE,
    show_colnames = FALSE,
    legend_breaks = distance_legend$breaks,
    legend_labels = distance_legend$labels,
    main = "Sample distance heatmap (Euclidean distance on blind VST)"
  )
  grDevices::dev.off()
}

# Marker panels are kept as simple named lookup tables so the gene sets are
# explicit in the script and easy to extend later.
build_sorting_marker_reference <- function() {
  tibble::tibble(
    panel_name = "sorting_markers",
    target_population = c("NeuN", "SOX10", "SOX2", "PU1"),
    requested_gene_name = c("RBFOX3", "SOX10", "SOX2", "SPI1"),
    marker_display_name = c("RBFOX3", "SOX10", "SOX2", "SPI1"),
    panel_plot_title = "Sorting marker gene expression across sorted populations"
  )
}

build_canonical_marker_reference <- function() {
  tibble::tibble(
    panel_name = "canonical_markers",
    target_population = c(
      "NeuN", "NeuN", "NeuN", "NeuN", "NeuN",
      "SOX2", "SOX2", "SOX2",
      "SOX10", "SOX10", "SOX10",
      "PU1", "PU1", "PU1",
      NA_character_, NA_character_, NA_character_
    ),
    requested_gene_name = c(
      "MAP2", "ENO2", "SNAP25", "TUBB3", "UCHL1",
      "GFAP", "SLC1A3", "ALDH1L1",
      "MBP", "MOG", "OPALIN",
      "TMEM119", "ITGAM", "SALL1",
      "PRNP", "GAL3ST1", "STX6"
    ),
    marker_display_name = c(
      "MAP2", "ENO2", "SNAP25", "TUBB3 (TUJ1)", "UCHL1 (PGP9.5)",
      "GFAP", "SLC1A3", "ALDH1L1",
      "MBP", "MOG", "OPALIN",
      "TMEM119", "ITGAM (CD11b)", "SALL1",
      "PRNP", "GAL3ST1", "STX6"
    ),
    target_population_display = c(
      rep(NA_character_, 17)
    ),
    panel_plot_title = "Canonical lineage marker expression across sorted populations"
  )
}

build_marker_gene_table <- function(gene_annotation, marker_reference) {
  if (!"target_population_display" %in% colnames(marker_reference)) {
    marker_reference$target_population_display <- NA_character_
  }

  marker_reference <- marker_reference |>
    dplyr::mutate(
      gene_name_upper = stringr::str_to_upper(.data$requested_gene_name),
      target_population_display = dplyr::if_else(
        is.na(.data$target_population_display),
        unname(exp383_population_display_labels[.data$target_population]),
        .data$target_population_display
      ),
      marker_order = dplyr::row_number()
    )

  annotation_lookup <- gene_annotation |>
    dplyr::mutate(gene_name_upper = stringr::str_to_upper(.data$gene_name)) |>
    dplyr::distinct(.data$gene_name_upper, .keep_all = TRUE)

  marker_reference |>
    dplyr::left_join(
      annotation_lookup,
      by = "gene_name_upper"
    ) |>
    dplyr::mutate(
      marker_label = dplyr::if_else(
        is.na(.data$target_population_display) | .data$target_population_display == "",
        .data$marker_display_name,
        paste0(.data$target_population_display, " | ", .data$marker_display_name)
      ),
      marker_found = !is.na(.data$gene_id)
    ) |>
    dplyr::select(
      panel_name,
      target_population,
      target_population_display,
      requested_gene_name,
      marker_display_name,
      matched_gene_name = gene_name,
      gene_id,
      marker_label,
      marker_order,
      panel_plot_title,
      marker_found
    )
}

build_marker_expression_table <- function(vsd, coldata, marker_table) {
  marker_found <- marker_table |>
    dplyr::filter(.data$marker_found)

  if (nrow(marker_found) == 0) {
    return(tibble::tibble())
  }

  marker_matrix <- SummarizedExperiment::assay(vsd)[marker_found$gene_id, , drop = FALSE]
  rownames(marker_matrix) <- marker_found$marker_label

  marker_long <- as.data.frame(marker_matrix) |>
    tibble::rownames_to_column("marker_label") |>
    tidyr::pivot_longer(
      cols = -marker_label,
      names_to = "sample",
      values_to = "vst_expression"
    ) |>
    dplyr::left_join(
      marker_found |>
        dplyr::select(
          panel_name,
          marker_label,
          marker_order,
          panel_plot_title,
          target_population,
          requested_gene_name,
          marker_display_name,
          matched_gene_name,
          gene_id
        ),
      by = "marker_label"
    ) |>
    dplyr::left_join(
      metadata_to_tibble(coldata),
      by = "sample"
    ) |>
    dplyr::mutate(
      population = exp383_population_factor(.data$population),
      target_population = exp383_population_factor(.data$target_population),
      marker_label = factor(
        .data$marker_label,
        levels = marker_table |>
          dplyr::arrange(.data$marker_order) |>
          dplyr::pull(.data$marker_label)
      )
    )

  marker_long
}

save_marker_expression_plot <- function(marker_expression, output_file, ncol = 2, width = 9, height = 7) {
  if (nrow(marker_expression) == 0) {
    return(invisible(NULL))
  }

  # Keep all facets on a common y-range so marker abundance differences are
  # directly comparable at a glance.
  global_limits <- range(marker_expression$vst_expression, na.rm = TRUE)
  plot_title <- unique(marker_expression$panel_plot_title)
  if (length(plot_title) != 1) {
    plot_title <- "Marker gene expression across sorted populations"
  }

  p <- ggplot2::ggplot(
    marker_expression,
    ggplot2::aes(x = .data$population, y = .data$vst_expression, colour = .data$population)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.25, width = 0.7) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.55, size = 1.2) +
    ggplot2::facet_wrap(~ marker_label, ncol = ncol) +
    exp383_scale_colour_population(name = "population") +
    exp383_scale_x_population() +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.05))
    ) +
    ggplot2::coord_cartesian(ylim = global_limits) +
    exp383_theme(base_size = 12) +
    ggplot2::labs(
      title = plot_title,
      x = "Sorted nuclei population",
      y = "Normalised expression"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  ggplot2::ggsave(output_file, plot = p, width = width, height = height)
}

save_marker_heatmap <- function(vsd, coldata, marker_table, output_file) {
  marker_found <- marker_table |>
    dplyr::filter(.data$marker_found)

  if (nrow(marker_found) == 0) {
    return(invisible(NULL))
  }

  sample_order <- coldata |>
    dplyr::mutate(sample = rownames(coldata)) |>
    dplyr::arrange(.data$population, .data$group_assignment, .data$sample) |>
    dplyr::pull(.data$sample)

  marker_matrix <- SummarizedExperiment::assay(vsd)[marker_found$gene_id, sample_order, drop = FALSE]
  rownames(marker_matrix) <- marker_found$marker_label

  # Row-scale the heatmap so that each marker's relative pattern across samples
  # is easier to see, even when genes differ in absolute abundance.
  marker_matrix_scaled <- t(scale(t(marker_matrix)))
  marker_limit <- max(abs(marker_matrix_scaled), na.rm = TRUE)
  marker_breaks <- seq(-marker_limit, marker_limit, length.out = 101)
  marker_legend_breaks <- c(-marker_limit, 0, marker_limit)
  marker_legend_labels <- format(round(marker_legend_breaks, 1), trim = TRUE, nsmall = 1)

  # As above, only annotate population. The purpose here is to check whether
  # sorted populations separate cleanly, not to visualise high-cardinality
  # experimental-group structure.
  annotation_df <- data.frame(
    population = exp383_population_display_factor(coldata[sample_order, "population"]),
    row.names = sample_order
  )

  grDevices::pdf(output_file, width = 12, height = 5)
  pheatmap::pheatmap(
    marker_matrix_scaled,
    annotation_col = annotation_df,
    annotation_colors = exp383_population_annotation_colors(),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_colnames = FALSE,
    color = grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
    breaks = marker_breaks,
    legend_breaks = marker_legend_breaks,
    legend_labels = marker_legend_labels,
    main = "FANS marker heatmap (row-scaled z-scores from blind VST)",
    name = "Z score"
  )
  grDevices::dev.off()
}

# -----------------
# Main script logic
# -----------------
args <- parse_cli_args()
stop_if_missing(args$input_rds, "Prepared input RDS")

analysis_input <- readRDS(args$input_rds)
required_objects <- c("txi_gene", "sample_metadata", "gene_annotation")
missing_objects <- setdiff(required_objects, names(analysis_input))
if (length(missing_objects) > 0) {
  stop(
    "Input RDS does not look like the output of 01_build_salmon_gene_inputs.R. ",
    "Missing objects: ", paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

coldata <- analysis_input$sample_metadata |>
  dplyr::mutate(
    population = exp383_population_factor(.data$population)
  ) |>
  as.data.frame()

rownames(coldata) <- coldata$sample

txi_gene <- analysis_input$txi_gene
sample_order <- rownames(coldata)
sample_idx <- match(sample_order, colnames(txi_gene$counts))
if (any(is.na(sample_idx))) {
  stop("Sample metadata and tximport counts are not aligned.", call. = FALSE)
}

txi_gene <- subset_txi(
  txi = txi_gene,
  keep_rows = rep(TRUE, nrow(txi_gene$counts)),
  keep_cols = sample_idx
)

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(args$output_dir, "plots")
table_dir <- file.path(args$output_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

keep_genes <- compute_grouped_filter(
  count_matrix = txi_gene$counts,
  grouping = coldata[[args$filter_group_var]],
  min_count = args$min_count,
  min_samples = args$min_samples
)

txi_filtered <- subset_txi(
  txi = txi_gene,
  keep_rows = keep_genes,
  keep_cols = rep(TRUE, ncol(txi_gene$counts))
)

filter_summary <- tibble::tibble(
  n_samples = nrow(coldata),
  n_genes_before_filter = nrow(txi_gene$counts),
  n_genes_after_filter = nrow(txi_filtered$counts),
  min_count = args$min_count,
  min_samples = args$min_samples,
  filter_group_var = args$filter_group_var,
  pca_ntop = args$pca_ntop
)
write_tsv(filter_summary, file.path(table_dir, "filter_summary.tsv"))
write_tsv(metadata_to_tibble(coldata), file.path(table_dir, "sample_metadata_used.tsv"))

dds <- DESeq2::DESeqDataSetFromTximport(
  txi = txi_filtered,
  colData = coldata,
  design = stats::as.formula("~ 1")
)
dds <- DESeq2::estimateSizeFactors(dds)
vsd <- DESeq2::vst(dds, blind = TRUE)

size_factors <- tibble::tibble(
  sample = colnames(dds),
  size_factor = DESeq2::sizeFactors(dds)
)
write_tsv(size_factors, file.path(table_dir, "size_factors.tsv"))

pca_table <- build_pca_coordinates(vsd, coldata, ntop = args$pca_ntop)
percent_variance <- attr(pca_table, "percent_variance")

pca_variance <- tibble::tibble(
  pc = paste0("PC", seq_along(percent_variance)),
  variance_percent = round(percent_variance, 4),
  cumulative_variance_percent = round(cumsum(percent_variance), 4)
)
write_tsv(pca_variance, file.path(table_dir, "pca_variance.tsv"))
write_tsv(pca_table, file.path(table_dir, "pca_coordinates.tsv"))

save_pca_plot(
  pca_table = pca_table,
  color_var = "population",
  percent_variance = percent_variance,
  output_file = file.path(plot_dir, "global_pca_by_population.pdf"),
  ntop = args$pca_ntop
)
save_pca_plot(
  pca_table = pca_table,
  color_var = "group_assignment",
  percent_variance = percent_variance,
  output_file = file.path(plot_dir, "global_pca_by_group_assignment.pdf"),
  ntop = args$pca_ntop
)
save_pca_plot(
  pca_table = pca_table,
  color_var = "inoculation_batch",
  percent_variance = percent_variance,
  output_file = file.path(plot_dir, "global_pca_by_inoculation_batch.pdf"),
  ntop = args$pca_ntop
)
save_pca_plot(
  pca_table = pca_table,
  color_var = "salmon_percent_mapped",
  percent_variance = percent_variance,
  output_file = file.path(plot_dir, "global_pca_by_salmon_percent_mapped.pdf"),
  ntop = args$pca_ntop
)
save_pca_plot(
  pca_table = pca_table,
  color_var = "rrna_alignment_percent",
  percent_variance = percent_variance,
  output_file = file.path(plot_dir, "global_pca_by_rrna_alignment_percent.pdf"),
  ntop = args$pca_ntop
)

save_sample_distance_heatmap(
  vsd = vsd,
  coldata = coldata,
  output_file = file.path(plot_dir, "global_sample_distance_heatmap.pdf")
)

sorting_marker_table <- build_marker_gene_table(
  gene_annotation = analysis_input$gene_annotation,
  marker_reference = build_sorting_marker_reference()
)
write_tsv(sorting_marker_table, file.path(table_dir, "marker_gene_lookup.tsv"))

sorting_marker_expression <- build_marker_expression_table(vsd, coldata, sorting_marker_table)
if (nrow(sorting_marker_expression) > 0) {
  write_tsv(sorting_marker_expression, file.path(table_dir, "marker_expression_vst.tsv"))
}

save_marker_expression_plot(
  marker_expression = sorting_marker_expression,
  output_file = file.path(plot_dir, "fans_marker_expression_by_population.pdf"),
  ncol = 2,
  width = 9,
  height = 7
)
save_marker_heatmap(
  vsd = vsd,
  coldata = coldata,
  marker_table = sorting_marker_table,
  output_file = file.path(plot_dir, "fans_marker_heatmap.pdf")
)

canonical_marker_table <- build_marker_gene_table(
  gene_annotation = analysis_input$gene_annotation,
  marker_reference = build_canonical_marker_reference()
)
write_tsv(canonical_marker_table, file.path(table_dir, "canonical_marker_gene_lookup.tsv"))

canonical_marker_expression <- build_marker_expression_table(vsd, coldata, canonical_marker_table)
if (nrow(canonical_marker_expression) > 0) {
  write_tsv(canonical_marker_expression, file.path(table_dir, "canonical_marker_expression_vst.tsv"))
}

save_marker_expression_plot(
  marker_expression = canonical_marker_expression,
  output_file = file.path(plot_dir, "canonical_marker_expression_by_population.pdf"),
  ncol = 3,
  width = 12,
  height = 13
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(args$output_dir, "sessionInfo.txt")
)

message("Finished global all-sample overview")
message("Samples used: ", nrow(coldata))
message("Genes retained after filtering: ", nrow(txi_filtered$counts))
