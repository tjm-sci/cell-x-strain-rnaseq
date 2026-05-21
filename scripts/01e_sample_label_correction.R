#!/usr/bin/env Rscript

# Global PCA indicates possible sample label mixup between a SOX2 and PU1 sample
# from 917442, consistent with the fact populations from any single animal had
# RNA extracted into adjacent tubes. Decision made to manually re-label these
# two samples.
#
# This script is intentionally narrow. It does not infer or search for other
# swaps. It applies the documented correction, preserves the original sample
# labels in metadata columns, and writes a corrected input object for downstream
# covariate screening and DEA.

suppressPackageStartupMessages({
  required_packages <- c("here")
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

suppressMessages(here::i_am("scripts/01e_sample_label_correction.R"))
source(here::here("scripts", "path_helpers.R"))

parse_cli_args <- function() {
  defaults <- list(
    input_rds = "results/dea/01_build_salmon_gene_inputs/exp383_salmon_gene_input.rds",
    output_dir = "results/dea/01e_sample_label_correction"
  )

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %% 2 != 0) {
    stop("Arguments must be provided as --key value pairs.", call. = FALSE)
  }

  parsed <- defaults
  if (length(args) > 0) {
    for (i in seq(1, length(args), by = 2)) {
      key <- sub("^--", "", args[[i]])
      value <- args[[i + 1]]

      if (!key %in% names(parsed)) {
        stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
      }

      parsed[[key]] <- value
    }
  }

  parsed
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
}

write_tsv <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

corrections <- data.frame(
  original_sample = c("917442-ME7-60-2-SOX2", "917442-ME7-60-2-PU1"),
  corrected_sample = c("917442-ME7-60-2-PU1", "917442-ME7-60-2-SOX2"),
  original_population = c("SOX2", "PU1"),
  corrected_population = c("PU1", "SOX2"),
  correction_reason = paste(
    "Manual relabel after 01b global PCA indicated the two adjacent-tube",
    "917442 SOX2 and PU1 samples aligned with the opposite labelled population."
  ),
  stringsAsFactors = FALSE
)

make_sample_map <- function(corrections) {
  sample_map <- stats::setNames(corrections$corrected_sample, corrections$original_sample)
  population_map <- stats::setNames(corrections$corrected_population, corrections$original_sample)
  list(sample = sample_map, population = population_map)
}

apply_sample_correction_to_table <- function(x, corrections, table_name, population_col = NULL) {
  if (!"sample" %in% colnames(x)) {
    stop(sprintf("%s does not contain a sample column.", table_name), call. = FALSE)
  }

  provenance_cols <- c("sample_label_original", "sample_label_correction_applied", "sample_label_correction_reason")
  if (!is.null(population_col)) {
    provenance_cols <- c(provenance_cols, "population_label_original")
  }

  existing_provenance <- intersect(provenance_cols, colnames(x))
  if (length(existing_provenance) > 0) {
    stop(
      sprintf("%s already contains correction provenance columns: %s", table_name, paste(existing_provenance, collapse = ", ")),
      call. = FALSE
    )
  }

  maps <- make_sample_map(corrections)
  correction_idx <- match(x$sample, corrections$original_sample)
  correction_applied <- !is.na(correction_idx)

  if (!all(corrections$original_sample %in% x$sample)) {
    missing_samples <- setdiff(corrections$original_sample, x$sample)
    stop(
      sprintf("%s is missing expected samples: %s", table_name, paste(missing_samples, collapse = ", ")),
      call. = FALSE
    )
  }

  x$sample_label_original <- x$sample
  x$sample_label_correction_applied <- correction_applied
  x$sample_label_correction_reason <- NA_character_
  x$sample_label_correction_reason[correction_applied] <- corrections$correction_reason[correction_idx[correction_applied]]

  if (!is.null(population_col)) {
    if (!population_col %in% colnames(x)) {
      stop(sprintf("%s does not contain %s.", table_name, population_col), call. = FALSE)
    }

    expected_population <- corrections$original_population[correction_idx[correction_applied]]
    observed_population <- as.character(x[[population_col]][correction_applied])
    if (!identical(observed_population, expected_population)) {
      stop(
        sprintf(
          "%s has unexpected original populations for corrected samples. Observed: %s. Expected: %s.",
          table_name,
          paste(observed_population, collapse = ", "),
          paste(expected_population, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    x$population_label_original <- as.character(x[[population_col]])
    x[[population_col]][correction_applied] <- unname(maps$population[x$sample[correction_applied]])
  }

  x$sample[correction_applied] <- unname(maps$sample[x$sample[correction_applied]])

  if (anyDuplicated(x$sample) > 0) {
    stop(sprintf("%s has duplicated sample labels after correction.", table_name), call. = FALSE)
  }

  x
}

apply_sample_correction_to_txi <- function(txi, corrections) {
  for (object_name in names(txi)) {
    object <- txi[[object_name]]
    if (is.null(dim(object)) || is.null(colnames(object))) {
      next
    }

    if (!all(corrections$original_sample %in% colnames(object))) {
      missing_samples <- setdiff(corrections$original_sample, colnames(object))
      stop(
        sprintf("txi_gene$%s is missing expected samples: %s", object_name, paste(missing_samples, collapse = ", ")),
        call. = FALSE
      )
    }

    corrected_colnames <- colnames(object)
    idx <- match(corrections$original_sample, corrected_colnames)
    corrected_colnames[idx] <- corrections$corrected_sample

    if (anyDuplicated(corrected_colnames) > 0) {
      stop(sprintf("txi_gene$%s has duplicated column names after correction.", object_name), call. = FALSE)
    }

    colnames(object) <- corrected_colnames
    txi[[object_name]] <- object
  }

  txi
}

args <- parse_cli_args()
input_rds <- resolve_project_path(args$input_rds, must_work = TRUE)
output_dir <- resolve_project_path(args$output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing(input_rds, "Input RDS")
analysis_input <- readRDS(input_rds)

required_objects <- c("txi_gene", "sample_metadata", "sample_manifest", "sample_qc_metrics")
missing_objects <- setdiff(required_objects, names(analysis_input))
if (length(missing_objects) > 0) {
  stop(sprintf("Input RDS is missing objects: %s", paste(missing_objects, collapse = ", ")), call. = FALSE)
}

analysis_input$txi_gene <- apply_sample_correction_to_txi(analysis_input$txi_gene, corrections)
analysis_input$sample_metadata <- apply_sample_correction_to_table(
  analysis_input$sample_metadata,
  corrections,
  table_name = "sample_metadata",
  population_col = "population"
)
analysis_input$sample_manifest <- apply_sample_correction_to_table(
  analysis_input$sample_manifest,
  corrections,
  table_name = "sample_manifest"
)
analysis_input$sample_qc_metrics <- apply_sample_correction_to_table(
  analysis_input$sample_qc_metrics,
  corrections,
  table_name = "sample_qc_metrics"
)

if (!identical(colnames(analysis_input$txi_gene$counts), analysis_input$sample_metadata$sample)) {
  stop("Corrected txi_gene column names do not match corrected sample_metadata$sample.", call. = FALSE)
}

analysis_input$settings$sample_label_correction <- list(
  applied = TRUE,
  correction_script = "scripts/01e_sample_label_correction.R",
  evidence_source = "results/dea/01b_global_sample_overview/tables/flagged_population_pca_samples.tsv",
  corrections = corrections
)

audit <- merge(
  corrections,
  analysis_input$sample_metadata[, c("sample", "sample_label_original", "fastq_1", "fastq_2")],
  by.x = "corrected_sample",
  by.y = "sample",
  all.x = TRUE,
  sort = FALSE
)
audit <- merge(
  audit,
  analysis_input$sample_manifest[, c("sample", "sample_label_original", "batch", "quant_sf")],
  by.x = "corrected_sample",
  by.y = "sample",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("_metadata", "_manifest")
)

output_rds <- file.path(output_dir, "exp383_salmon_gene_input_label_corrected.rds")
saveRDS(analysis_input, output_rds)
write_tsv(audit, file.path(output_dir, "sample_label_corrections_applied.tsv"))
write_tsv(analysis_input$sample_metadata, file.path(output_dir, "sample_metadata_label_corrected.tsv"))
write_tsv(corrections, file.path(output_dir, "manual_sample_label_corrections.tsv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message("Wrote corrected input object: ", output_rds)
