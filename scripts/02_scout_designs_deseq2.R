#!/usr/bin/env Rscript

# This script is the design-scouting stage for DEA.
#
# It expects the unified object produced by 01_build_salmon_gene_inputs.R and is
# intended to answer questions such as:
#
# - Should we use ~ group_assignment ?
# - Does adding inoculation_batch help or just consume degrees of freedom?
# - Do sample_mass_mg or incubation_time_hrs look like nuisance covariates worth
#   including?
#
# The script is intentionally separate from the final DEA script so we can try
# multiple candidate designs, inspect QC plots, and settle on a final model
# before producing inferential result tables.
#
# Output layout is explicitly design-specific:
#   results/dea/02_design_scout/<population>/<design_id>/plots
#   results/dea/02_design_scout/<population>/<design_id>/tables
#
# Shared, population-level artefacts that do not depend on design are written to:
#   results/dea/02_design_scout/<population>/shared

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

# Parse either a one-line design or a TSV with multiple candidate designs.
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

# The metadata were cleaned in script 01, but here we make sure any remaining
# character columns are converted to factors before they enter a model formula.
prepare_coldata_for_modeling <- function(coldata) {
  is_character <- vapply(coldata, is.character, logical(1))
  coldata[is_character] <- lapply(coldata[is_character], factor)
  coldata
}

# Filtering is done once per population and is intentionally independent of the
# specific design formula being tested. The rule is:
#
# Keep a gene if it has at least `min_count` estimated counts in at least
# `min_samples` libraries within at least one level of `filter_group_var`.
#
# For this study, the default grouping variable is `group_assignment`, which is a
# sensible way to avoid keeping genes that never rise above background within any
# biological condition.
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
  percent_var <- round(100 * attr(pca, "percentVar"), 2)
  pca$sample <- rownames(pca)

  # Add the remaining metadata columns so any covariate can be plotted later.
  merge(pca, cbind(sample = rownames(coldata), as.data.frame(coldata)), by = "sample", sort = FALSE)
}

save_pca_plot <- function(pca_table, color_var, title_text, output_file) {
  if (!color_var %in% colnames(pca_table)) {
    return(invisible(NULL))
  }

  plot_data <- pca_table
  plot_data[[color_var]] <- plot_data[[color_var]]

  percent_var <- attr(pca_table, "percent_var")
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC1, y = PC2, color = .data[[color_var]])) +
    ggplot2::geom_point(size = 2.8, alpha = 0.9) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = title_text,
      x = sprintf("PC1 (%.2f%% variance)", percent_var[[1]]),
      y = sprintf("PC2 (%.2f%% variance)", percent_var[[2]]),
      color = color_var
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(output_file, plot = p, width = 7, height = 5.5)
}

save_sample_distance_heatmap <- function(vsd, annotation_df, output_file) {
  sample_dist <- dist(t(SummarizedExperiment::assay(vsd)))
  sample_dist_matrix <- as.matrix(sample_dist)

  gr_devices <- grDevices::pdf
  gr_devices(output_file, width = 9, height = 8)
  pheatmap::pheatmap(
    sample_dist_matrix,
    annotation_col = annotation_df,
    annotation_row = annotation_df,
    main = "Sample-to-sample distances"
  )
  grDevices::dev.off()
}

save_dispersion_plot <- function(dds, output_file) {
  grDevices::pdf(output_file, width = 7, height = 5.5)
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
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      title = "DESeq2 size factors",
      x = "Sample",
      y = "Size factor"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(output_file, plot = p, width = 7, height = 10)
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

if (!args$population %in% as.character(unique(coldata_all$population))) {
  stop(
    sprintf("Population '%s' not found in prepared metadata.", args$population),
    call. = FALSE
  )
}

# Subset to the selected sorted population before filtering and modeling.
keep_samples <- coldata_all$population == args$population
coldata_pop <- coldata_all[keep_samples, , drop = FALSE]
rownames(coldata_pop) <- coldata_pop$sample

txi_pop <- subset_txi(
  txi = txi_gene,
  keep_rows = rep(TRUE, nrow(txi_gene$counts)),
  keep_cols = keep_samples
)

# Choose a sensible default min_samples from the smallest group size in the
# selected population. This makes the default filter align with the actual
# replication available for that population.
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

# Evaluate every candidate design separately, but always on the same population
# subset and filtered gene set.
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

  # Record the exact design text for traceability.
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

  # Build a DESeqDataSet directly from the tximport output so gene-length offsets
  # are carried forward in the standard DESeq2-aware way.
  dds <- DESeq2::DESeqDataSetFromTximport(
    txi = txi_pop_filtered,
    colData = coldata_pop,
    design = design_formula_obj
  )

  # Run the standard DESeq2 fitting sequence. We are not extracting contrasts in
  # this script; the point here is to estimate size factors and dispersions, then
  # inspect how the design behaves.
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  vsd <- DESeq2::vst(dds, blind = FALSE)

  normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
  size_factor_table <- data.frame(
    sample = colnames(dds),
    size_factor = DESeq2::sizeFactors(dds),
    group_assignment = coldata_pop$group_assignment,
    stringsAsFactors = FALSE
  )

  pca_table <- build_pca_table(vsd, coldata_pop)
  attr(pca_table, "percent_var") <- round(100 * attr(DESeq2::plotPCA(vsd, intgroup = "group_assignment", returnData = TRUE), "percentVar"), 2)

  # Store the design matrix because it is often the quickest way to sanity-check
  # whether a candidate formula is encoded as intended.
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

  # Core QC plots that should always be present for a design scout.
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

  # Produce PCA plots for the design terms themselves, plus a small standard set
  # of high-value covariates for quick interpretation.
  variables_to_plot <- unique(c(
    "group_assignment",
    all.vars(design_formula_obj),
    "inoculum",
    "dpi_factor",
    "inoculation_batch",
    "sample_mass_mg",
    "incubation_time_hrs",
    "date_nuc_prep_factor"
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
