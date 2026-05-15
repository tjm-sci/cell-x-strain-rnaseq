#!/usr/bin/env Rscript

# This script quantifies how similar each EXP383 sample is to purified
# Zhang et al. 2014 cortex reference populations. The contamination-focused
# summaries still emphasize astrocytes, OPCs, and oligodendrocytes, but the
# joint reference panel also includes neurons and microglia as extra context:
#   1. purified astrocytes
#   2. purified OPCs
#   3. purified oligodendrocytes
#   4. purified neurons
#   5. purified microglia
#
# It is a follow-on exploratory script for the SOX2 contamination question. The
# goal is not to estimate literal contamination fractions, but to provide a
# numeric and visual description of lineage similarity using the shared marker
# panel derived in 01c.
#
# Important distinction:
# - the numeric centroid-similarity scores use the matched Zhang marker panel
#   from 01c, because that panel was built specifically to ask a contamination
#   question about astrocyte / OPC / oligodendrocyte identity
# - the joint PCA does NOT use that marker panel; instead it uses a broader
#   top-variable-gene space derived from all Zhang reference genes so the
#   visualization is less biased by the glial contamination screen
#
# The workflow is split into three stages:
# 1. build the Zhang reference centroids in the 01c marker space
# 2. standardize each dataset gene-wise and calculate centroid similarities
#    from that marker space
# 3. build a separate 500-gene Zhang-derived PCA space for visualization,
#    showing reference centroids as large square points over the EXP383 sample
#    cloud

suppressPackageStartupMessages({
  required_packages <- c(
    "DESeq2",
    "dplyr",
    "ggplot2",
    "ggrepel",
    "tibble",
    "tidyr",
    "here"
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

# Reuse the project-wide theme and population colours so this exploratory output
# looks consistent with the rest of the repo.
suppressMessages(here::i_am("scripts/01d_sox2_reference_similarity.R"))
source(here::here("scripts", "path_helpers.R"))
source(here::here("scripts", "plot_style.R"))

reference_display_levels <- c(
  "Purified astrocytes",
  "Purified OPCs",
  "Purified oligodendrocytes",
  "Purified neurons",
  "Purified microglia"
)

### COMMAND-LINE ARGUMENTS ####################################################

# Keep the CLI small. This script is project-specific and builds directly on the
# cached outputs from 01b and 01c.
parse_cli_args <- function() {
  defaults <- list(
    overview_rds = "results/dea/01b_global_sample_overview/cache/global_overview_input.rds",
    input_rds = "results/dea/01_build_salmon_gene_inputs/exp383_salmon_gene_input.rds",
    marker_csv = "results/dea/01c_sox2_contamination_eda/tables/zhang_cell_type_enriched_markers.csv",
    reference_avg_csv = "results/dea/01c_sox2_contamination_eda/tables/zhang_reference_cell_type_averages.csv",
    output_dir = "results/dea/01d_sox2_reference_similarity",
    pca_ntop = "5000"
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

  parsed$pca_ntop <- as.integer(parsed$pca_ntop)

  parsed
}

### SMALL HELPERS #############################################################

# Stop early when a required file is missing.
stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
}

# Write tidy outputs as CSV for straightforward review outside R.
write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, quote = TRUE)
}

# Convert metadata to a tibble while preserving the canonical sample identifier.
metadata_to_tibble <- function(coldata) {
  metadata_tbl <- tibble::as_tibble(coldata)

  if (!"sample" %in% colnames(metadata_tbl)) {
    metadata_tbl <- tibble::rownames_to_column(metadata_tbl, "sample")
  }

  metadata_tbl
}

# Standardize each gene within a matrix so sample-to-reference comparison is
# driven by relative expression pattern rather than absolute scale. The same
# operation is applied separately to the EXP383 sample matrix and the Zhang
# reference centroids.
row_zscore_matrix <- function(mat) {
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds <- apply(mat, 1, stats::sd, na.rm = TRUE)
  row_sds[is.na(row_sds) | row_sds == 0] <- 1

  scaled <- (mat - row_means) / row_sds
  rownames(scaled) <- rownames(mat)
  colnames(scaled) <- colnames(mat)
  scaled
}

### LOAD THE CACHED 01B OVERVIEW INPUT ########################################

# 01d uses the same global filter and blind VST that already underpin the 01b
# overview, so we load that cache rather than rebuilding the expression object.
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

# The broader joint PCA needs the full EXP383 gene annotation so Zhang gene
# symbols can be matched outside the narrower marker panel used for the
# contamination-oriented similarity scores.
read_analysis_input <- function(input_rds) {
  analysis_input <- readRDS(input_rds)
  required_objects <- c("gene_annotation")
  missing_objects <- setdiff(required_objects, names(analysis_input))

  if (length(missing_objects) > 0) {
    stop(
      "Input RDS does not look like the output of 01_build_salmon_gene_inputs.R. ",
      "Missing objects: ", paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  analysis_input
}

### 1. BUILD THE SHARED REFERENCE CENTROIDS ###################################

# Load the full matched marker set from 01c. This table defines the feature
# space used for the numeric centroid-similarity scores. Unlike the plotted 01c
# boxplot, we use all matched markers that passed the Zhang enrichment rule, not
# just the top six per cell type.
read_marker_panel <- function(marker_csv) {
  utils::read.csv(marker_csv) |>
    tibble::as_tibble() |>
    dplyr::filter(.data$exp383_gene_found) |>
    dplyr::select(
      "target_reference_cell_type",
      "target_cell_type_display",
      "gene_symbol",
      "exp383_gene_id",
      "exp383_gene_name",
      "fold_enrichment_vs_other_mean",
      "fold_enrichment_vs_other_max"
    ) |>
    dplyr::distinct()
}

# Load the Zhang cell-type averages created in 01c and keep the purified
# reference classes used in this similarity analysis. These same averages feed
# two different downstream paths:
# - the full matched-marker centroid matrix used for numeric similarity scores
# - the broader all-gene matrix later used to derive the joint PCA feature set
read_reference_averages <- function(reference_avg_csv) {
  utils::read.csv(reference_avg_csv) |>
    tibble::as_tibble() |>
    dplyr::filter(
      .data$reference_cell_type %in% c(
        "astrocytes",
        "opc",
        "myelinating_oligodendrocyte",
        "neurons",
        "microglia_macrophage"
      )
    )
}

# Build one gene-by-reference centroid matrix in the matched 01c marker space.
# This object is used only for the numeric centroid-similarity calculations.
# The joint PCA deliberately uses a different, broader feature set later.
build_reference_centroid_matrix <- function(marker_panel, reference_averages) {
  marker_panel |>
    dplyr::select("gene_symbol", "exp383_gene_id", "exp383_gene_name") |>
    dplyr::distinct() |>
    dplyr::inner_join(
      reference_averages,
      by = "gene_symbol"
    ) |>
    dplyr::select("gene_symbol", "exp383_gene_id", "exp383_gene_name", "reference_cell_type", "mean_reference_expression") |>
    tidyr::pivot_wider(
      names_from = "reference_cell_type",
      values_from = "mean_reference_expression"
    ) |>
    dplyr::rename(
      "Purified astrocytes" = "astrocytes",
      "Purified OPCs" = "opc",
      "Purified oligodendrocytes" = "myelinating_oligodendrocyte",
      "Purified neurons" = "neurons",
      "Purified microglia" = "microglia_macrophage"
    )
}

# Keep only the 01c markers that truly exist in the cached 01b VST matrix. This
# makes the marker-based similarity space explicit and avoids carrying forward a
# small number of markers that passed 01c matching but were dropped by the
# global 01b count filter.
align_reference_centroids_to_vsd <- function(reference_centroids, vsd) {
  available_gene_ids <- rownames(SummarizedExperiment::assay(vsd))

  reference_centroids |>
    dplyr::filter(.data$exp383_gene_id %in% available_gene_ids)
}

### 2. STANDARDIZE AND CALCULATE CENTROID SIMILARITY ##########################

# Build the EXP383 sample matrix on the same 01c-derived marker genes used for
# the reference centroids. This matrix underpins the numeric sample-to-centroid
# similarity scores, not the later joint PCA.
build_exp383_sample_matrix <- function(vsd, reference_centroids) {
  assay_matrix <- SummarizedExperiment::assay(vsd)
  match_idx <- match(reference_centroids$exp383_gene_id, rownames(assay_matrix))

  sample_matrix <- assay_matrix[match_idx, , drop = FALSE]
  rownames(sample_matrix) <- reference_centroids$gene_symbol
  sample_matrix
}

# Calculate sample-to-centroid Pearson correlations after gene-wise
# standardization in the 01c marker space. These are relative similarity scores,
# not contamination percentages, and they are intentionally separate from the
# broader feature space used for the joint PCA.
calculate_reference_similarity <- function(sample_matrix, reference_matrix, coldata) {
  sample_scaled <- row_zscore_matrix(sample_matrix)
  reference_scaled <- row_zscore_matrix(reference_matrix)

  similarity_matrix <- stats::cor(sample_scaled, reference_scaled, method = "pearson")

  similarity_long <- as.data.frame(similarity_matrix) |>
    tibble::rownames_to_column("sample") |>
    tidyr::pivot_longer(
      cols = -sample,
      names_to = "reference_cell_type",
      values_to = "pearson_similarity"
    ) |>
    dplyr::left_join(
      metadata_to_tibble(coldata),
      by = "sample"
    ) |>
    dplyr::mutate(
      population = exp383_population_factor(.data$population),
      reference_cell_type = factor(
        .data$reference_cell_type,
        levels = reference_display_levels
      )
    )

  similarity_wide <- similarity_long |>
    dplyr::select("sample", "population", "inoculum", "dpi", "group_assignment", "reference_cell_type", "pearson_similarity") |>
    tidyr::pivot_wider(
      names_from = "reference_cell_type",
      values_from = "pearson_similarity"
    ) |>
    dplyr::mutate(
      opc_minus_astro = .data$`Purified OPCs` - .data$`Purified astrocytes`,
      oligodendrocyte_minus_astro = .data$`Purified oligodendrocytes` - .data$`Purified astrocytes`,
      oligolineage_max_minus_astro = pmax(.data$`Purified OPCs`, .data$`Purified oligodendrocytes`) - .data$`Purified astrocytes`
    )

  list(
    similarity_long = similarity_long,
    similarity_wide = similarity_wide,
    sample_scaled = sample_scaled,
    reference_scaled = reference_scaled
  )
}

# Summarize the similarity structure with a population-wise boxplot.
save_similarity_boxplot <- function(similarity_long, output_file) {
  p <- ggplot2::ggplot(
    similarity_long,
    ggplot2::aes(x = .data$population, y = .data$pearson_similarity, colour = .data$population)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.25, width = 0.7) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.55, size = 1.2) +
    ggplot2::facet_wrap(~ reference_cell_type, ncol = 3) +
    exp383_scale_colour_population(name = "population") +
    exp383_scale_x_population() +
    ggplot2::coord_cartesian(ylim = c(-1, 1)) +
    exp383_theme(base_size = 12) +
    ggplot2::labs(
      title = "Zhang reference centroid similarity by sorted population",
      subtitle = "Pearson correlation after gene-wise z-scoring within dataset",
      x = "Sorted nuclei population",
      y = "Reference centroid similarity (Pearson r)"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  exp383_save_ggplot(output_file, plot = p, width = 12, height = 6)
}

### 3. JOINT PCA VISUALIZATION ################################################

# For the joint PCA, switch to a broader and less contamination-biased feature
# space than the one used for centroid scoring:
# - start from all genes in the Zhang reference matrix
# - log-transform Zhang reference expression
# - identify the top `ntop` most variable genes across the purified reference
#   cell types
# - identify the top `ntop` most variable genes across the cached EXP383 VST
#   matrix
# - match genes across datasets and keep the intersection only
#
# This means the PCA geometry is not driven by the narrower 01c marker panel.
build_joint_pca_gene_space <- function(reference_averages, gene_annotation, vsd, ntop) {
  reference_wide <- reference_averages |>
    dplyr::mutate(reference_log_expression = log2(.data$mean_reference_expression + 1)) |>
    dplyr::select("gene_symbol", "reference_cell_type", "reference_log_expression") |>
    tidyr::pivot_wider(
      names_from = "reference_cell_type",
      values_from = "reference_log_expression"
    ) |>
    dplyr::rename(
      "Purified astrocytes" = "astrocytes",
      "Purified OPCs" = "opc",
      "Purified oligodendrocytes" = "myelinating_oligodendrocyte",
      "Purified neurons" = "neurons",
      "Purified microglia" = "microglia_macrophage"
    )

  reference_matrix <- reference_wide |>
    dplyr::select(dplyr::all_of(reference_display_levels)) |>
    as.matrix()
  rownames(reference_matrix) <- reference_wide$gene_symbol

  reference_variance <- apply(reference_matrix, 1, stats::var, na.rm = TRUE)
  reference_ranked <- reference_wide |>
    dplyr::mutate(
      reference_variance = reference_variance
    ) |>
    dplyr::arrange(dplyr::desc(.data$reference_variance), .data$gene_symbol)

  annotation_lookup <- gene_annotation |>
    dplyr::mutate(gene_name_upper = stringr::str_to_upper(.data$gene_name)) |>
    dplyr::distinct(.data$gene_name_upper, .keep_all = TRUE)

  assay_matrix <- SummarizedExperiment::assay(vsd)
  available_gene_ids <- rownames(assay_matrix)
  exp383_variance <- apply(assay_matrix, 1, stats::var, na.rm = TRUE)

  exp383_top_gene_table <- tibble::tibble(
    exp383_gene_id = rownames(assay_matrix),
    exp383_variance = exp383_variance
  ) |>
    dplyr::arrange(dplyr::desc(.data$exp383_variance), .data$exp383_gene_id) |>
    dplyr::slice_head(n = ntop) |>
    dplyr::mutate(exp383_rank = dplyr::row_number())

  pca_gene_table <- reference_ranked |>
    dplyr::mutate(
      reference_rank = dplyr::row_number(),
      gene_name_upper = stringr::str_to_upper(.data$gene_symbol)
    ) |>
    dplyr::left_join(
      annotation_lookup,
      by = "gene_name_upper"
    ) |>
    dplyr::rename(
      exp383_gene_id = "gene_id",
      exp383_gene_name = "gene_name"
    ) |>
    dplyr::filter(.data$exp383_gene_id %in% available_gene_ids) |>
    dplyr::slice_head(n = ntop) |>
    dplyr::inner_join(
      exp383_top_gene_table,
      by = "exp383_gene_id"
    ) |>
    dplyr::mutate(
      combined_rank = .data$reference_rank + .data$exp383_rank
    ) |>
    dplyr::arrange(.data$combined_rank, dplyr::desc(.data$reference_variance), dplyr::desc(.data$exp383_variance), .data$gene_symbol)

  if (nrow(pca_gene_table) == 0) {
    stop(
      "The intersection of the top variable Zhang genes and top variable EXP383 genes is empty.",
      call. = FALSE
    )
  }

  reference_pca_matrix <- reference_matrix[pca_gene_table$gene_symbol, , drop = FALSE]

  match_idx <- match(pca_gene_table$exp383_gene_id, rownames(assay_matrix))
  sample_pca_matrix <- assay_matrix[match_idx, , drop = FALSE]
  rownames(sample_pca_matrix) <- pca_gene_table$gene_symbol

  list(
    pca_gene_table = pca_gene_table,
    sample_pca_matrix = sample_pca_matrix,
    reference_pca_matrix = reference_pca_matrix
  )
}

# Run a joint PCA on the standardized broad-gene matrix built above. This is a
# visualization-only embedding and should be interpreted separately from the
# marker-based centroid-similarity scores.
run_joint_pca <- function(sample_scaled, reference_scaled, coldata) {
  combined_matrix <- cbind(sample_scaled, reference_scaled)
  pca <- stats::prcomp(t(combined_matrix), center = TRUE, scale. = FALSE)
  percent_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

  combined_scores <- as.data.frame(pca$x) |>
    tibble::rownames_to_column("observation") |>
    dplyr::mutate(
      source = dplyr::if_else(
        .data$observation %in% colnames(reference_scaled),
        "Zhang reference centroid",
        "EXP383 sample"
      )
    )

  sample_scores <- combined_scores |>
    dplyr::filter(.data$source == "EXP383 sample") |>
    dplyr::rename(sample = "observation") |>
    dplyr::left_join(
      metadata_to_tibble(coldata),
      by = "sample"
    ) |>
    dplyr::mutate(
      population = exp383_population_factor(.data$population)
    )

  reference_scores <- combined_scores |>
    dplyr::filter(.data$source == "Zhang reference centroid") |>
    dplyr::rename(reference_cell_type = "observation") |>
    dplyr::mutate(
      reference_cell_type = factor(
        .data$reference_cell_type,
        levels = reference_display_levels
      )
    )

  list(
    sample_scores = sample_scores,
    reference_scores = reference_scores,
    percent_var = percent_var
  )
}

# Plot the joint PCA with the reference centroids drawn last as large square
# markers so they stand out clearly from the EXP383 sample cloud.
save_joint_pca_plot <- function(sample_scores, reference_scores, percent_var, output_file, ntop) {
  reference_palette <- c(
    "Purified astrocytes" = unname(exp383_population_palette[["SOX2"]]),
    "Purified OPCs" = "#b35806",
    "Purified oligodendrocytes" = unname(exp383_population_palette[["SOX10"]]),
    "Purified neurons" = unname(exp383_population_palette[["NeuN"]]),
    "Purified microglia" = unname(exp383_population_palette[["PU1"]])
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = sample_scores,
      mapping = ggplot2::aes(x = .data$PC1, y = .data$PC2, colour = .data$population),
      size = 2.2,
      alpha = 0.75
    ) +
    ggplot2::geom_point(
      data = reference_scores,
      mapping = ggplot2::aes(x = .data$PC1, y = .data$PC2, fill = .data$reference_cell_type),
      shape = 22,
      size = 5.5,
      stroke = 0.8,
      colour = "black"
    ) +
    ggrepel::geom_label_repel(
      data = reference_scores,
      mapping = ggplot2::aes(x = .data$PC1, y = .data$PC2, label = .data$reference_cell_type),
      seed = 383,
      size = 3,
      fill = "white",
      colour = "black",
      label.size = 0.2,
      box.padding = 0.35,
      point.padding = 0.3,
      segment.color = "grey35",
      segment.size = 0.3,
      max.overlaps = Inf,
      min.segment.length = 0
    ) +
    exp383_scale_colour_population(name = "EXP383 population") +
    ggplot2::scale_fill_manual(values = reference_palette, name = "Zhang reference") +
    exp383_theme(base_size = 12) +
    ggplot2::labs(
      title = "Joint PCA of EXP383 samples and Zhang reference centroids",
      subtitle = sprintf(
        "Intersection of the top %d most variable genes in EXP383 (blind VST) and Zhang (log2 reference means), then gene-wise z-scoring within dataset",
        ntop
      ),
      x = sprintf("PC1 (%.2f%% variance)", percent_var[[1]]),
      y = sprintf("PC2 (%.2f%% variance)", percent_var[[2]])
    )

  exp383_save_ggplot(output_file, plot = p, width = 8, height = 6)
}

### MAIN SCRIPT LOGIC #########################################################

args <- parse_cli_args()
args$overview_rds <- resolve_project_path(args$overview_rds)
args$input_rds <- resolve_project_path(args$input_rds)
args$marker_csv <- resolve_project_path(args$marker_csv)
args$reference_avg_csv <- resolve_project_path(args$reference_avg_csv)
args$output_dir <- resolve_project_path(args$output_dir)
stop_if_missing(args$overview_rds, "01b cached overview RDS")
stop_if_missing(args$input_rds, "01_build input RDS")
stop_if_missing(args$marker_csv, "01c marker CSV")
stop_if_missing(args$reference_avg_csv, "01c Zhang reference average CSV")

overview_input <- read_overview_input(args$overview_rds)
analysis_input <- read_analysis_input(args$input_rds)
marker_panel <- read_marker_panel(args$marker_csv)
reference_averages <- read_reference_averages(args$reference_avg_csv)

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(args$output_dir, "plots")
table_dir <- file.path(args$output_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

reference_centroids <- build_reference_centroid_matrix(
  marker_panel = marker_panel,
  reference_averages = reference_averages
)
reference_centroids <- align_reference_centroids_to_vsd(
  reference_centroids = reference_centroids,
  vsd = overview_input$vsd
)
write_csv(reference_centroids, file.path(table_dir, "zhang_reference_centroids_shared_marker_space.csv"))

reference_matrix <- reference_centroids |>
  dplyr::select(dplyr::all_of(reference_display_levels)) |>
  as.matrix()
rownames(reference_matrix) <- reference_centroids$gene_symbol

sample_matrix <- build_exp383_sample_matrix(
  vsd = overview_input$vsd,
  reference_centroids = reference_centroids
)

similarity_results <- calculate_reference_similarity(
  sample_matrix = sample_matrix,
  reference_matrix = reference_matrix,
  coldata = overview_input$coldata
)
write_csv(
  similarity_results$similarity_long,
  file.path(table_dir, "reference_similarity_scores_long.csv")
)
write_csv(
  similarity_results$similarity_wide,
  file.path(table_dir, "reference_similarity_scores_wide.csv")
)
write_csv(
  dplyr::filter(similarity_results$similarity_wide, .data$population == "SOX2") |>
    dplyr::arrange(dplyr::desc(.data$oligolineage_max_minus_astro)),
  file.path(table_dir, "sox2_reference_similarity_scores.csv")
)

save_similarity_boxplot(
  similarity_long = similarity_results$similarity_long,
  output_file = file.path(plot_dir, "reference_centroid_similarity_by_population.png")
)

joint_pca_gene_space <- build_joint_pca_gene_space(
  reference_averages = reference_averages,
  gene_annotation = analysis_input$gene_annotation,
  vsd = overview_input$vsd,
  ntop = args$pca_ntop
)
write_csv(
  joint_pca_gene_space$pca_gene_table,
  file.path(table_dir, "joint_pca_gene_panel.csv")
)

joint_pca <- run_joint_pca(
  sample_scaled = row_zscore_matrix(joint_pca_gene_space$sample_pca_matrix),
  reference_scaled = row_zscore_matrix(joint_pca_gene_space$reference_pca_matrix),
  coldata = overview_input$coldata
)
write_csv(
  joint_pca$sample_scores,
  file.path(table_dir, "joint_pca_sample_coordinates.csv")
)
write_csv(
  joint_pca$reference_scores,
  file.path(table_dir, "joint_pca_reference_coordinates.csv")
)
write_csv(
  tibble::tibble(
    pc = paste0("PC", seq_along(joint_pca$percent_var)),
    variance_percent = joint_pca$percent_var,
    cumulative_variance_percent = cumsum(joint_pca$percent_var)
  ),
  file.path(table_dir, "joint_pca_variance.csv")
)

save_joint_pca_plot(
  sample_scores = joint_pca$sample_scores,
  reference_scores = joint_pca$reference_scores,
  percent_var = joint_pca$percent_var,
  output_file = file.path(plot_dir, "joint_reference_similarity_pca.png"),
  ntop = args$pca_ntop
)

run_summary <- tibble::tibble(
  n_samples = nrow(overview_input$coldata),
  n_genes_after_filter = nrow(overview_input$txi_filtered$counts),
  n_markers_in_01c_panel = nrow(marker_panel),
  n_shared_markers = nrow(reference_centroids),
  n_reference_cell_types = ncol(reference_matrix),
  joint_pca_ntop = args$pca_ntop,
  joint_pca_shared_genes = nrow(joint_pca_gene_space$pca_gene_table)
)
write_csv(run_summary, file.path(table_dir, "run_summary.csv"))

message("Finished SOX2 reference similarity analysis")
message("Samples used: ", nrow(overview_input$coldata))
message("Shared markers used: ", nrow(reference_centroids))
