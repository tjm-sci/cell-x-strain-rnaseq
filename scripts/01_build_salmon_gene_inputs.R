#!/usr/bin/env Rscript

# This script builds a single, analysis-ready gene-level input object for all
# 286 libraries in EXP383.
#
# The nf-core/rnaseq run was executed in four separate population batches, so
# the Salmon outputs currently live in four different result directories. Before
# we can do any downstream DEA, we want one place where:
#
# 1. every sample has its Salmon quant.sf file recorded
# 2. metadata are cleaned and aligned to those quant files
# 3. sample-level QC metrics from nf-core/MultiQC are attached once
# 4. tximport has been run once in a reproducible way
# 5. the resulting object is saved for later scripts
#
# This script deliberately does NOT perform DEA, differential filtering, or
# model fitting. Those steps are left to later scripts so that design formulas
# can be changed without rebuilding the input object from scratch.

suppressPackageStartupMessages({
  # We only need tximport here. All other file handling uses base R.
  if (!requireNamespace("tximport", quietly = TRUE)) {
    stop(
      "Package 'tximport' is required but is not installed. ",
      "Install it in your R environment before running this script.",
      call. = FALSE
    )
  }
})

# -----------------------------
# Command-line argument parsing
# -----------------------------
# To keep this script portable, argument parsing is done with a very small base-R
# helper instead of requiring optparse/argparse.
parse_cli_args <- function() {
  defaults <- list(
    results_root = "/media/tmurphy/4TB_HDD/exp383/nfcore_rnaseq",
    metadata_csv = "metadata/exp383_nfcore_rnaseq_samplesheet_template.csv",
    output_dir = "results/dea/01_build_salmon_gene_inputs"
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

# -----------------
# Utility functions
# -----------------
stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
}

write_tsv <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

# The sample sheet already contains most covariates we need, but we make a few
# types explicit here so later DEA scripts do not need to guess them.
clean_metadata <- function(metadata) {
  required_columns <- c(
    "sample", "mouse_n", "population", "inoculum", "dpi", "group_assignment",
    "inoculation_batch", "sample_mass_mg", "date_nuc_prep",
    "incubation_time_hrs"
  )

  missing_columns <- setdiff(required_columns, colnames(metadata))
  if (length(missing_columns) > 0) {
    stop(
      sprintf(
        "Metadata file is missing required columns: %s",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  metadata$sample <- as.character(metadata$sample)
  metadata$mouse_n <- as.character(metadata$mouse_n)
  metadata$population <- factor(
    metadata$population,
    levels = c("NeuN", "SOX10", "SOX2", "PU1")
  )
  metadata$inoculum <- factor(metadata$inoculum)
  metadata$group_assignment <- factor(metadata$group_assignment)
  metadata$inoculation_batch <- factor(metadata$inoculation_batch)

  # Keep both numeric and factor versions of dpi because different downstream
  # models may want different representations.
  metadata$dpi <- as.integer(metadata$dpi)
  metadata$dpi_factor <- factor(metadata$dpi, levels = sort(unique(metadata$dpi)))

  # Technical covariates that are naturally numeric are parsed now so later
  # scripts can use them directly in model formulas.
  metadata$sample_mass_mg <- suppressWarnings(as.numeric(metadata$sample_mass_mg))
  metadata$incubation_time_hrs <- suppressWarnings(as.numeric(metadata$incubation_time_hrs))

  # Dates are kept in three forms:
  # - the original string column
  # - a factor for categorical modeling
  # - an integer day offset for numeric trend modeling
  prep_dates <- as.Date(metadata$date_nuc_prep)
  metadata$date_nuc_prep_factor <- factor(metadata$date_nuc_prep)
  if (all(is.na(prep_dates))) {
    metadata$date_nuc_prep_days <- NA_integer_
  } else {
    metadata$date_nuc_prep_days <- as.integer(prep_dates - min(prep_dates, na.rm = TRUE))
  }

  # Provide centered/scaled numeric covariates once here so later design files
  # can use stable, model-friendly column names instead of repeating ad hoc
  # scaling inside every analysis script.
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

  metadata$sample_mass_mg_z <- scale_numeric_covariate(metadata$sample_mass_mg)
  metadata$incubation_time_hrs_z <- scale_numeric_covariate(metadata$incubation_time_hrs)
  metadata$date_nuc_prep_days_z <- scale_numeric_covariate(metadata$date_nuc_prep_days)

  metadata
}

# Salmon quant files live under:
#   <results_root>/<batch>/salmon/<sample>/quant.sf
# We discover them rather than hard-coding per-batch paths.
discover_quant_files <- function(results_root) {
  quant_files <- list.files(
    path = results_root,
    pattern = "quant\\.sf$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(quant_files) == 0) {
    stop("No Salmon quant.sf files were found under results_root.", call. = FALSE)
  }

  manifest <- data.frame(
    sample = basename(dirname(quant_files)),
    batch = basename(dirname(dirname(dirname(quant_files)))),
    quant_sf = normalizePath(quant_files, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )

  duplicate_samples <- unique(manifest$sample[duplicated(manifest$sample)])
  if (length(duplicate_samples) > 0) {
    stop(
      sprintf(
        "Duplicate sample IDs were discovered in quant files: %s",
        paste(duplicate_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  manifest[order(manifest$sample), , drop = FALSE]
}

# nf-core writes a compact, tabular summary of per-sample QC metrics into each
# batch MultiQC directory. These are exactly the kind of technical covariates we
# want to screen before final DEA, so we import them once here rather than
# rebuilding them later in every downstream script.
discover_sample_qc_metrics <- function(results_root) {
  general_stats_files <- list.files(
    path = results_root,
    pattern = "multiqc_general_stats\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(general_stats_files) == 0) {
    stop("No multiqc_general_stats.txt files were found under results_root.", call. = FALSE)
  }

  selected_columns <- c(
    "Sample",
    "salmon-percent_mapped",
    "salmon-num_mapped",
    "fastqc_raw-percent_duplicates",
    "fastqc_raw-percent_gc",
    "fastqc_raw-total_sequences",
    "fastp-pct_duplication",
    "fastp-after_filtering_q30_rate",
    "fastp-filtering_result_passed_filter_reads",
    "fastp-after_filtering_gc_content",
    "fastp-pct_surviving",
    "fastp-pct_adapter",
    "bowtie2_rrna_removal-overall_alignment_rate",
    "fastqc_filtered-percent_duplicates",
    "fastqc_filtered-percent_gc",
    "fastqc_filtered-total_sequences"
  )

  rename_map <- c(
    "Sample" = "sample",
    "salmon-percent_mapped" = "salmon_percent_mapped",
    "salmon-num_mapped" = "salmon_num_mapped_millions",
    "fastqc_raw-percent_duplicates" = "raw_percent_duplicates",
    "fastqc_raw-percent_gc" = "raw_percent_gc",
    "fastqc_raw-total_sequences" = "raw_total_sequences_millions",
    "fastp-pct_duplication" = "fastp_percent_duplication",
    "fastp-after_filtering_q30_rate" = "fastp_q30_rate_after_filtering",
    "fastp-filtering_result_passed_filter_reads" = "fastp_passed_filter_reads_millions",
    "fastp-after_filtering_gc_content" = "fastp_percent_gc_after_filtering",
    "fastp-pct_surviving" = "fastp_percent_surviving",
    "fastp-pct_adapter" = "fastp_percent_adapter",
    "bowtie2_rrna_removal-overall_alignment_rate" = "rrna_alignment_percent",
    "fastqc_filtered-percent_duplicates" = "filtered_percent_duplicates",
    "fastqc_filtered-percent_gc" = "filtered_percent_gc",
    "fastqc_filtered-total_sequences" = "filtered_total_sequences_millions"
  )

  qc_tables <- lapply(general_stats_files, function(path) {
    qc <- read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

    missing_columns <- setdiff(selected_columns, colnames(qc))
    if (length(missing_columns) > 0) {
      stop(
        sprintf(
          "MultiQC general stats file is missing expected columns:\n- %s\nMissing: %s",
          path,
          paste(missing_columns, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    qc <- qc[, selected_columns, drop = FALSE]

    # MultiQC includes per-read companion rows such as "Sample Read 1".
    # Downstream DEA covariates should only use the sample-level summary rows.
    qc <- qc[!grepl(" Read [12]$", qc$Sample), , drop = FALSE]
    colnames(qc) <- unname(rename_map[colnames(qc)])

    numeric_columns <- setdiff(colnames(qc), "sample")
    qc[numeric_columns] <- lapply(qc[numeric_columns], function(x) suppressWarnings(as.numeric(x)))

    qc
  })

  qc_metrics <- do.call(rbind, qc_tables)
  duplicate_samples <- unique(qc_metrics$sample[duplicated(qc_metrics$sample)])
  if (length(duplicate_samples) > 0) {
    stop(
      sprintf(
        "Duplicate sample IDs were discovered across MultiQC general stats files: %s",
        paste(duplicate_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  qc_metrics[order(qc_metrics$sample), , drop = FALSE]
}

# The tx2gene mapping should be the same across all four batch outputs because
# they all used the same Ensembl reference. We verify that assumption explicitly.
load_and_check_tx2gene <- function(results_root) {
  tx2gene_files <- list.files(
    path = results_root,
    pattern = "salmon\\.merged\\.tx2gene\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(tx2gene_files) == 0) {
    stop("No salmon.merged.tx2gene.tsv files were found.", call. = FALSE)
  }

  tx2gene_tables <- lapply(tx2gene_files, function(path) {
    read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  })

  reference_table <- tx2gene_tables[[1]]
  for (i in seq_along(tx2gene_tables)[-1]) {
    if (!identical(reference_table, tx2gene_tables[[i]])) {
      stop(
        sprintf(
          "tx2gene mismatch detected between batch files:\n- %s\n- %s",
          tx2gene_files[[1]], tx2gene_files[[i]]
        ),
        call. = FALSE
      )
    }
  }

  reference_table
}

# -----------------
# Main script logic
# -----------------
args <- parse_cli_args()
results_root <- normalizePath(args$results_root, winslash = "/", mustWork = FALSE)
metadata_csv <- normalizePath(args$metadata_csv, winslash = "/", mustWork = FALSE)
output_dir <- args$output_dir

stop_if_missing(results_root, "Results root")
stop_if_missing(metadata_csv, "Metadata CSV")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading metadata: ", metadata_csv)
metadata <- read.csv(metadata_csv, stringsAsFactors = FALSE, check.names = FALSE)
metadata <- clean_metadata(metadata)

message("Discovering Salmon quant files under: ", results_root)
manifest <- discover_quant_files(results_root)

message("Discovering sample-level MultiQC metrics under: ", results_root)
qc_metrics <- discover_sample_qc_metrics(results_root)

expected_samples <- metadata$sample
observed_samples <- manifest$sample

missing_quant <- setdiff(expected_samples, observed_samples)
extra_quant <- setdiff(observed_samples, expected_samples)

if (length(missing_quant) > 0) {
  stop(
    sprintf(
      "The following samples are present in metadata but missing quant.sf files: %s",
      paste(missing_quant, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (length(extra_quant) > 0) {
  stop(
    sprintf(
      "The following quant.sf samples are not present in metadata: %s",
      paste(extra_quant, collapse = ", ")
    ),
    call. = FALSE
  )
}

missing_qc <- setdiff(expected_samples, qc_metrics$sample)
extra_qc <- setdiff(qc_metrics$sample, expected_samples)

if (length(missing_qc) > 0) {
  stop(
    sprintf(
      "The following samples are present in metadata but missing MultiQC sample metrics: %s",
      paste(missing_qc, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (length(extra_qc) > 0) {
  stop(
    sprintf(
      "The following MultiQC sample metrics are not present in metadata: %s",
      paste(extra_qc, collapse = ", ")
    ),
    call. = FALSE
  )
}

# Reorder the manifest to the metadata sample order so every downstream matrix
# has the same column order as the metadata table.
manifest <- manifest[match(metadata$sample, manifest$sample), , drop = FALSE]
stopifnot(identical(metadata$sample, manifest$sample))

qc_metrics <- qc_metrics[match(metadata$sample, qc_metrics$sample), , drop = FALSE]
stopifnot(identical(metadata$sample, qc_metrics$sample))

# Fold sequencing/QC covariates directly into the sample metadata so later DEA
# scripts only need one metadata table.
metadata <- cbind(
  metadata,
  qc_metrics[, setdiff(colnames(qc_metrics), "sample"), drop = FALSE]
)

message("Loading and checking tx2gene tables")
tx2gene <- load_and_check_tx2gene(results_root)

message("Running tximport across all samples")
quant_files_named <- manifest$quant_sf
names(quant_files_named) <- manifest$sample

txi_gene <- tximport::tximport(
  files = quant_files_named,
  type = "salmon",
  tx2gene = tx2gene[, c("transcript_id", "gene_id")],
  countsFromAbundance = "no",
  # We only need the core abundance/count/length matrices for DEA.
  # Any inferential replicate metadata shipped by Salmon can be dropped here.
  dropInfReps = TRUE,
  ignoreTxVersion = FALSE
)

# Gene name annotations are helpful later when writing result tables.
gene_annotation <- unique(tx2gene[, c("gene_id", "gene_name")])
gene_annotation <- gene_annotation[order(gene_annotation$gene_id), , drop = FALSE]

# This object is the main handoff to later DEA scripts.
analysis_input <- list(
  txi_gene = txi_gene,
  sample_metadata = metadata,
  sample_manifest = manifest,
  sample_qc_metrics = qc_metrics,
  tx2gene = tx2gene,
  gene_annotation = gene_annotation,
  settings = list(
    countsFromAbundance = "no",
    created_at = as.character(Sys.time()),
    results_root = results_root,
    metadata_csv = metadata_csv
  )
)

message("Writing output files to: ", output_dir)
saveRDS(analysis_input, file = file.path(output_dir, "exp383_salmon_gene_input.rds"))
write_tsv(metadata, file.path(output_dir, "sample_metadata_cleaned.tsv"))
write_tsv(manifest, file.path(output_dir, "sample_to_quant_manifest.tsv"))
write_tsv(qc_metrics, file.path(output_dir, "sample_qc_metrics.tsv"))
write_tsv(tx2gene, file.path(output_dir, "tx2gene.tsv"))
write_tsv(gene_annotation, file.path(output_dir, "gene_annotation.tsv"))

population_counts <- as.data.frame(table(metadata$population), stringsAsFactors = FALSE)
colnames(population_counts) <- c("population", "n_samples")
write_tsv(population_counts, file.path(output_dir, "population_sample_counts.tsv"))

writeLines(capture.output(sessionInfo()), con = file.path(output_dir, "sessionInfo.txt"))

message("Finished building unified Salmon gene-level input object")
message("Samples imported: ", nrow(metadata))
message("Genes imported: ", nrow(txi_gene$counts))
message("Batches merged: ", paste(sort(unique(manifest$batch)), collapse = ", "))
