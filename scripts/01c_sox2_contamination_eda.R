#!/usr/bin/env Rscript

# This script uses the Zhang et al. 2014 purified mouse cortex RNA-seq resource
# to derive reference marker sets for three glial cell classes relevant to the
# SOX2 contamination question:
#   1. purified astrocytes
#   2. purified OPCs
#   3. purified oligodendrocytes
#
# The workflow is deliberately simple and specific to this project:
# - average the duplicate Zhang replicates per reference cell type
# - extract a gene symbol from the Zhang gene label
# - define enriched markers using an explicit fold-enrichment rule
# - match those markers into the EXP383 RNA-seq object
# - plot their VST expression across the four sorted populations
#
# This is exploratory contamination-focused analysis only. It does not perform
# any differential testing.

suppressPackageStartupMessages({
  required_packages <- c(
    "dplyr",
    "ggplot2",
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

### SHARED PLOTTING HELPERS ###################################################

# Reuse the project-wide colour palette and ggplot theme so these exploratory
# figures sit naturally beside the 01b global overview outputs.
get_script_dir <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 0) {
    return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  }

  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE))
}

source(file.path(get_script_dir(), "plot_style.R"))

### COMMAND-LINE ARGUMENTS ####################################################

# Keep the CLI small. This script is meant to be reused inside this repo, not
# as a general-purpose marker-discovery tool.
parse_cli_args <- function() {
  defaults <- list(
    overview_rds = "results/dea/01b_global_sample_overview/cache/global_overview_input.rds",
    input_rds = "results/dea/01_build_salmon_gene_inputs/exp383_salmon_gene_input.rds",
    zhang_csv = "resources/zhang-et-al-2014-mouse-rnaseq.csv",
    output_dir = "results/dea/01c_sox2_contamination_eda",
    min_target_expression = "20",
    fold_change_cutoff = "5",
    max_other_fold_cutoff = "2",
    plot_top_n = "6"
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
  parsed$min_target_expression <- as.numeric(parsed$min_target_expression)
  parsed$fold_change_cutoff <- as.numeric(parsed$fold_change_cutoff)
  parsed$max_other_fold_cutoff <- as.numeric(parsed$max_other_fold_cutoff)
  parsed$plot_top_n <- as.integer(parsed$plot_top_n)

  parsed
}

### SMALL HELPERS #############################################################

# Stop early when a required input is missing.
stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
}

# Write tidy tabular outputs as CSV so they are easy to open in spreadsheet
# software during exploratory review.
write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, quote = TRUE)
}

# Convert the sample metadata to a tibble without losing the canonical sample
# identifier.
metadata_to_tibble <- function(coldata) {
  metadata_tbl <- tibble::as_tibble(coldata)

  if (!"sample" %in% colnames(metadata_tbl)) {
    metadata_tbl <- tibble::rownames_to_column(metadata_tbl, "sample")
  }

  metadata_tbl
}
### LOAD THE CACHED GLOBAL OVERVIEW INPUT #####################################

# 01c asks exactly the same question on the same transformed expression space as
# 01b, so it is cleaner to reuse the cached overview handoff rather than
# rebuilding the global filter and VST matrix.
read_overview_input <- function(overview_rds) {
  overview_input <- readRDS(overview_rds)
  required_objects <- c("vsd", "coldata", "txi_filtered", "filter_summary", "settings")
  missing_objects <- setdiff(required_objects, names(overview_input))

  if (length(missing_objects) > 0) {
    stop(
      "Overview RDS does not look like the output cached by 01b_global_sample_overview.R. ",
      "Missing objects: ", paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  overview_input
}

### WRANGLE THE ZHANG REFERENCE MATRIX ########################################

# The downloaded Zhang file stores duplicate purified-cell replicates in wide
# format. We convert it to a tidy format and add a clean gene-symbol column.
read_zhang_reference <- function(zhang_csv) {
  zhang_tbl <- utils::read.csv(zhang_csv, check.names = FALSE) |>
    tibble::as_tibble()

  zhang_tbl |>
    dplyr::mutate(
      gene_symbol = stringr::str_remove(.data$gene_id, " - Mus musculus$")
    ) |>
    tidyr::pivot_longer(
      cols = -c("gene_id", "id", "gene_symbol"),
      names_to = "sample_name",
      values_to = "reference_expression"
    ) |>
    dplyr::mutate(
      reference_cell_type = stringr::str_remove(.data$sample_name, "_[12]$")
    )
}

# Average the duplicate purified-cell replicates so each Zhang reference cell
# type contributes one representative expression profile per gene.
average_zhang_replicates <- function(zhang_long) {
  zhang_long |>
    dplyr::group_by(.data$gene_id, .data$gene_symbol, .data$reference_cell_type) |>
    dplyr::summarise(
      mean_reference_expression = mean(.data$reference_expression, na.rm = TRUE),
      .groups = "drop"
    )
}

# The contamination question focuses on three reference classes. The Zhang
# dataset also contains neurons, microglia/macrophages, endothelial cells, and
# newly formed oligodendrocytes; these remain in the background set when we
# compute enrichment.
build_target_class_reference <- function() {
  tibble::tibble(
    target_reference_cell_type = c(
      "astrocytes",
      "opc",
      "myelinating_oligodendrocyte"
    ),
    target_cell_type_display = c(
      "Purified astrocytes",
      "Purified OPCs",
      "Purified oligodendrocytes"
    )
  )
}

# Define a simple, explicit enrichment rule.
#
# A gene is called a marker for a target reference cell type when:
# - its mean expression in the target Zhang cell type is at least 20
# - it is enriched at least 5-fold over the mean of all other Zhang reference
#   cell types
# - and it is still enriched at least 2-fold over the single highest competing
#   reference cell type, i.e. the second-closest cell-type profile
#
# This keeps the rule readable while avoiding genes that only look specific
# because most non-target cell types are close to zero.
find_zhang_enriched_markers <- function(
  zhang_average,
  min_target_expression,
  fold_change_cutoff,
  max_other_fold_cutoff
) {
  zhang_wide <- zhang_average |>
    dplyr::select("gene_id", "gene_symbol", "reference_cell_type", "mean_reference_expression") |>
    tidyr::pivot_wider(
      names_from = "reference_cell_type",
      values_from = "mean_reference_expression"
    )

  target_classes <- build_target_class_reference()
  all_cell_type_columns <- setdiff(colnames(zhang_wide), c("gene_id", "gene_symbol"))

  marker_tables <- lapply(seq_len(nrow(target_classes)), function(i) {
    target_type <- target_classes$target_reference_cell_type[[i]]
    target_display <- target_classes$target_cell_type_display[[i]]
    other_types <- setdiff(all_cell_type_columns, target_type)

    target_expression <- zhang_wide[[target_type]]
    other_expression_matrix <- as.matrix(zhang_wide[, other_types, drop = FALSE])
    other_mean_expression <- rowMeans(other_expression_matrix, na.rm = TRUE)
    other_max_expression <- apply(other_expression_matrix, 1, max, na.rm = TRUE)

    zhang_wide |>
      dplyr::mutate(
        target_mean_expression = target_expression,
        other_mean_expression = other_mean_expression,
        other_max_expression = other_max_expression,
        fold_enrichment_vs_other_mean = (.data$target_mean_expression + 1) / (.data$other_mean_expression + 1),
        fold_enrichment_vs_other_max = (.data$target_mean_expression + 1) / (.data$other_max_expression + 1),
        target_reference_cell_type = target_type,
        target_cell_type_display = target_display
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        passes_marker_rule = (
          .data$target_mean_expression >= min_target_expression &
            .data$fold_enrichment_vs_other_mean >= fold_change_cutoff &
            .data$fold_enrichment_vs_other_max >= max_other_fold_cutoff
        )
      ) |>
      dplyr::filter(.data$passes_marker_rule) |>
      dplyr::arrange(dplyr::desc(.data$fold_enrichment_vs_other_mean), dplyr::desc(.data$target_mean_expression)) |>
      dplyr::mutate(
        rank_within_target = dplyr::row_number()
      )
  })

  dplyr::bind_rows(marker_tables)
}

### MATCH ZHANG MARKERS INTO THE EXP383 MATRIX ################################

# Match the Zhang-derived marker genes into the EXP383 gene annotation so we can
# extract their VST expression for plotting.
match_markers_to_exp383 <- function(marker_table, gene_annotation) {
  annotation_lookup <- gene_annotation |>
    dplyr::mutate(gene_name_upper = stringr::str_to_upper(.data$gene_name)) |>
    dplyr::distinct(.data$gene_name_upper, .keep_all = TRUE)

  marker_table |>
    dplyr::mutate(
      gene_name_upper = stringr::str_to_upper(.data$gene_symbol)
    ) |>
    dplyr::left_join(
      annotation_lookup,
      by = "gene_name_upper"
    ) |>
    dplyr::mutate(
      exp383_gene_found = !is.na(.data$gene_id.y)
    ) |>
    dplyr::rename_with(
      ~ c("zhang_gene_id", "exp383_gene_id", "exp383_gene_name"),
      .cols = c("gene_id.x", "gene_id.y", "gene_name")
    ) |>
    dplyr::select(
      "target_reference_cell_type",
      "target_cell_type_display",
      "gene_symbol",
      "zhang_gene_id",
      "target_mean_expression",
      "other_mean_expression",
      "other_max_expression",
      "fold_enrichment_vs_other_mean",
      "fold_enrichment_vs_other_max",
      "rank_within_target",
      "exp383_gene_id",
      "exp383_gene_name",
      "exp383_gene_found"
    )
}

# Build the compact plot panel using the strongest Zhang markers that are also
# present in the EXP383 expression object.
select_plot_markers <- function(marker_lookup, plot_top_n) {
  marker_lookup |>
    dplyr::filter(.data$exp383_gene_found) |>
    dplyr::group_by(.data$target_reference_cell_type) |>
    dplyr::arrange(
      dplyr::desc(.data$fold_enrichment_vs_other_mean),
      dplyr::desc(.data$target_mean_expression),
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = plot_top_n) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      marker_label = paste0(.data$target_cell_type_display, " | ", .data$gene_symbol),
      panel_plot_title = "Zhang et al. cell type enriched marker expression by population",
      marker_order = dplyr::row_number()
    )
}

### BUILD THE BOXPLOT INPUT TABLE #############################################

# Convert the selected marker set into a long VST expression table so we can
# make faceted population boxplots in the same style as script 01b.
build_marker_expression_table <- function(vsd, coldata, plot_markers) {
  if (nrow(plot_markers) == 0) {
    return(tibble::tibble())
  }

  marker_matrix <- SummarizedExperiment::assay(vsd)[plot_markers$exp383_gene_id, , drop = FALSE]
  rownames(marker_matrix) <- plot_markers$marker_label

  as.data.frame(marker_matrix) |>
    tibble::rownames_to_column("marker_label") |>
    tidyr::pivot_longer(
      cols = -marker_label,
      names_to = "sample",
      values_to = "vst_expression"
    ) |>
    dplyr::left_join(
      plot_markers |>
        dplyr::select(
          "target_reference_cell_type",
          "target_cell_type_display",
          "gene_symbol",
          "exp383_gene_id",
          "exp383_gene_name",
          "marker_label",
          "marker_order",
          "panel_plot_title"
        ),
      by = "marker_label"
    ) |>
    dplyr::left_join(
      metadata_to_tibble(coldata),
      by = "sample"
    ) |>
    dplyr::mutate(
      population = exp383_population_factor(.data$population),
      marker_label = factor(
        .data$marker_label,
        levels = plot_markers |>
          dplyr::arrange(.data$marker_order) |>
          dplyr::pull("marker_label")
      )
    )
}

### SAVE THE FACETED BOXPLOT ##################################################

# Use one fixed y-axis range across all facets so the astrocyte, OPC, and
# oligodendrocyte marker panels can be compared at a glance.
#
# The plotting table already contains the top six matched markers per purified
# Zhang reference cell type, ranked by enrichment over the mean of the other
# reference cell types and then by target-class abundance.
save_marker_expression_plot <- function(marker_expression, output_file) {
  if (nrow(marker_expression) == 0) {
    return(invisible(NULL))
  }

  global_limits <- range(marker_expression$vst_expression, na.rm = TRUE)
  plot_title <- unique(marker_expression$panel_plot_title)
  if (length(plot_title) != 1) {
    plot_title <- "Zhang et al. cell type enriched marker expression by population"
  }

  p <- ggplot2::ggplot(
    marker_expression,
    ggplot2::aes(x = .data$population, y = .data$vst_expression, colour = .data$population)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.25, width = 0.7) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.55, size = 1.2) +
    ggplot2::facet_wrap(~ marker_label, ncol = 3) +
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
      y = "Normalised expression (VST)"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  exp383_save_ggplot(output_file, plot = p, width = 12, height = 10)
}

### MAIN SCRIPT LOGIC #########################################################

args <- parse_cli_args()
stop_if_missing(args$overview_rds, "01b cached overview RDS")
stop_if_missing(args$zhang_csv, "Zhang reference CSV")

analysis_input <- readRDS(args$input_rds)
required_objects <- c("gene_annotation")
missing_objects <- setdiff(required_objects, names(analysis_input))
if (length(missing_objects) > 0) {
  stop(
    "Input RDS does not look like the output of 01_build_salmon_gene_inputs.R. ",
    "Missing objects: ", paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

overview_input <- read_overview_input(args$overview_rds)

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(args$output_dir, "plots")
table_dir <- file.path(args$output_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

zhang_long <- read_zhang_reference(args$zhang_csv)
zhang_average <- average_zhang_replicates(zhang_long)
write_csv(
  zhang_average,
  file.path(table_dir, "zhang_reference_cell_type_averages.csv")
)

zhang_markers <- find_zhang_enriched_markers(
  zhang_average = zhang_average,
  min_target_expression = args$min_target_expression,
  fold_change_cutoff = args$fold_change_cutoff,
  max_other_fold_cutoff = args$max_other_fold_cutoff
)

marker_lookup <- match_markers_to_exp383(
  marker_table = zhang_markers,
  gene_annotation = analysis_input$gene_annotation
)
write_csv(
  marker_lookup,
  file.path(table_dir, "zhang_cell_type_enriched_markers.csv")
)

plot_markers <- select_plot_markers(
  marker_lookup = marker_lookup,
  plot_top_n = args$plot_top_n
)
write_csv(
  plot_markers,
  file.path(table_dir, "zhang_cell_type_enriched_markers_for_plot.csv")
)

marker_expression <- build_marker_expression_table(
  vsd = overview_input$vsd,
  coldata = overview_input$coldata,
  plot_markers = plot_markers
)
write_csv(
  marker_expression,
  file.path(table_dir, "zhang_marker_expression_vst.csv")
)

save_marker_expression_plot(
  marker_expression = marker_expression,
  output_file = file.path(plot_dir, "zhang_cell_type_enriched_marker_expression_by_population.png")
)

filter_summary <- tibble::tibble(
  n_samples = nrow(overview_input$coldata),
  n_genes_after_filter = nrow(overview_input$txi_filtered$counts),
  zhang_marker_count = nrow(marker_lookup),
  plot_marker_count = nrow(plot_markers),
  overview_cache = args$overview_rds,
  min_target_expression = args$min_target_expression,
  fold_change_cutoff = args$fold_change_cutoff,
  max_other_fold_cutoff = args$max_other_fold_cutoff,
  plot_top_n = args$plot_top_n
)
write_csv(
  filter_summary,
  file.path(table_dir, "run_summary.csv")
)

message("Finished SOX2 contamination EDA")
message("Samples used: ", nrow(overview_input$coldata))
message("Genes retained after filtering: ", nrow(overview_input$txi_filtered$counts))
message("Zhang markers retained: ", nrow(marker_lookup))
message("Markers plotted: ", nrow(plot_markers))
