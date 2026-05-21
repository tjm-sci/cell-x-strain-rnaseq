#!/usr/bin/env Rscript

# This script builds the single output object used by downstream DEA scripts.
#
# The nf-core/rnaseq run was executed in four separate population batches, so
# the Salmon outputs, tx2gene tables, and MultiQC summaries currently live in
# four separate result directories. Before any downstream modelling, we want one
# stable object that:
#
# 1. records the `quant.sf` path for every sample
# 2. cleans and types the project metadata once so types are the same and correct
#    When we split the downstream analysis by sorted population.
# 3. joins one row of sample-level QC metrics per sample
# 4. runs tximport once in a reproducible way
# 5. saves the resulting object for all later scripts
#
# This script does not perform DEA or gene filtering. It just prepares the
# shared starting point.

suppressPackageStartupMessages({
  required_packages <- c("tximport", "dplyr", "tibble", "here")
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

suppressMessages(here::i_am("scripts/01_build_salmon_gene_inputs.R"))
source(here::here("scripts", "path_helpers.R"))

### COMMAND-LINE ARGUMENTS ####################################################


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

### SMALL HELPERS #############################################################

# Stop early when an expected file or directory is missing. This keeps the main
# workflow readable and keeps failures close to their real cause.
stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
}

load_excluded_mouse_ids <- function(repo_root) {
  exclusion_flags_file <- file.path(
    dirname(repo_root),
    "exp383_mouse_metadata",
    "output",
    "exp383_mouse_exclusion_flags.csv"
  )

  if (!file.exists(exclusion_flags_file)) {
    stop(
      "Project-level mouse exclusion flags not found: ",
      exclusion_flags_file,
      ". Run exp383_mouse_metadata/scripts/05_check_non_experimental_culls.R first.",
      call. = FALSE
    )
  }

  flags <- utils::read.csv(exclusion_flags_file, stringsAsFactors = FALSE, check.names = FALSE)
  required_columns <- c("animal_id", "exclude_from_downstream", "exclusion_applies_to")
  missing_columns <- setdiff(required_columns, colnames(flags))
  if (length(missing_columns) > 0) {
    stop(
      sprintf("Mouse exclusion flags missing columns: %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  flags |>
    dplyr::mutate(
      exclude_from_downstream = as.character(.data$exclude_from_downstream) %in% c("TRUE", "true", "1")
    ) |>
    dplyr::filter(
      .data$exclude_from_downstream,
      .data$exclusion_applies_to %in% c("all_exp383_repos", basename(repo_root))
    ) |>
    dplyr::pull(.data$animal_id) |>
    as.character() |>
    unique()
}

sample_mouse_ids <- function(sample_ids) {
  sub("-.*$", "", as.character(sample_ids))
}

# Write small tabular outputs in a consistent TSV format for later inspection.
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

# Standardize numeric covariates once here so downstream design formulas can use
# stable `_z` columns instead of repeating ad hoc scaling later.
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

### METADATA PREPARATION ######################################################

# The sample sheet already carries the core biological and technical metadata.
# Here we make the intended data types explicit so later DEA scripts do not need
# to guess whether a column should be numeric, factor, or both.
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

  metadata_tbl <- tibble::as_tibble(metadata)
  prep_dates <- as.Date(metadata_tbl$date_nuc_prep)

  date_nuc_prep_days <- if (all(is.na(prep_dates))) {
    rep(NA_integer_, nrow(metadata_tbl))
  } else {
    as.integer(prep_dates - min(prep_dates, na.rm = TRUE))
  }

  metadata_tbl |>
    dplyr::mutate(
      sample = as.character(.data$sample),
      mouse_n = as.character(.data$mouse_n),
      population = factor(.data$population, levels = c("NeuN", "SOX10", "SOX2", "PU1")),
      inoculum = factor(.data$inoculum),
      group_assignment = factor(.data$group_assignment),
      inoculation_batch = factor(.data$inoculation_batch),
      dpi = as.integer(.data$dpi),
      dpi_factor = factor(.data$dpi, levels = sort(unique(.data$dpi))),
      sample_mass_mg = suppressWarnings(as.numeric(.data$sample_mass_mg)),
      incubation_time_hrs = suppressWarnings(as.numeric(.data$incubation_time_hrs)),
      date_nuc_prep_factor = factor(.data$date_nuc_prep),
      date_nuc_prep_days = date_nuc_prep_days,
      sample_mass_mg_z = scale_numeric_covariate(.data$sample_mass_mg),
      incubation_time_hrs_z = scale_numeric_covariate(.data$incubation_time_hrs),
      date_nuc_prep_days_z = scale_numeric_covariate(.data$date_nuc_prep_days)
    )
}

### DISCOVER SALMON OUTPUTS ###################################################

# Salmon quant files live under:
#   <results_root>/<batch>/salmon/<sample>/quant.sf
# We discover them dynamically so this script remains aligned to the actual
# result tree on disk.
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

  manifest <- tibble::tibble(
    sample = basename(dirname(quant_files)),
    batch = basename(dirname(dirname(dirname(quant_files)))),
    quant_sf = normalizePath(quant_files, winslash = "/", mustWork = TRUE)
  ) |>
    dplyr::arrange(.data$sample)

  duplicate_samples <- manifest |>
    dplyr::count(.data$sample) |>
    dplyr::filter(.data$n > 1) |>
    dplyr::pull(.data$sample)

  if (length(duplicate_samples) > 0) {
    stop(
      sprintf(
        "Duplicate sample IDs were discovered in quant files: %s",
        paste(duplicate_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  manifest
}

### DISCOVER SAMPLE-LEVEL QC ##################################################

# nf-core writes a compact per-sample summary into each batch's
# `multiqc_general_stats.txt`. We import the sample-level rows once here and
# carry them forward into DEA as technical covariates.
#
# Important detail:
# - `fastp_percent_adapter` is the adapter burden detected and trimmed by
#   fastp from the input reads. It is not a post-filter residual adapter
#   percentage. MultiQC does not expose a comparable post-trim adapter column in
#   this summary table.
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
    qc <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

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

    qc <- qc |>
      tibble::as_tibble() |>
      dplyr::select(dplyr::all_of(selected_columns)) |>
      dplyr::filter(!grepl(" Read [12]$", .data$Sample)) |>
      stats::setNames(unname(rename_map[selected_columns]))

    numeric_columns <- setdiff(colnames(qc), "sample")
    qc |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(numeric_columns),
          ~ suppressWarnings(as.numeric(.x))
        )
      )
  })

  qc_metrics <- dplyr::bind_rows(qc_tables) |>
    dplyr::arrange(.data$sample)

  duplicate_samples <- qc_metrics |>
    dplyr::count(.data$sample) |>
    dplyr::filter(.data$n > 1) |>
    dplyr::pull(.data$sample)

  if (length(duplicate_samples) > 0) {
    stop(
      sprintf(
        "Duplicate sample IDs were discovered across MultiQC general stats files: %s",
        paste(duplicate_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  qc_metrics
}

### LOAD AND VERIFY TX2GENE ###################################################

# All four nf-core batches were run against the same Ensembl reference, so
# their tx2gene tables should be identical. We verify that explicitly here.
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
    utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE) |>
      tibble::as_tibble()
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

### MAIN WORKFLOW #############################################################

args <- parse_cli_args()
repo_root <- here::here()
results_root <- resolve_project_path(args$results_root)
metadata_csv <- resolve_project_path(args$metadata_csv)
output_dir <- resolve_project_path(args$output_dir)

stop_if_missing(results_root, "Results root")
stop_if_missing(metadata_csv, "Metadata CSV")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

excluded_mouse_ids <- load_excluded_mouse_ids(repo_root)

message("Reading metadata: ", metadata_csv)
metadata <- utils::read.csv(metadata_csv, stringsAsFactors = FALSE, check.names = FALSE) |>
  clean_metadata()

# Remove mice culled for non-experimental reasons.
# See ../exp383_mouse_metadata/output/exp383_mouse_exclusion_flags.csv.
metadata_row_count <- nrow(metadata)
metadata <- metadata |>
  dplyr::filter(!.data$mouse_n %in% excluded_mouse_ids)
message("Excluded RNA-seq metadata rows for non-experimental cull reasons: ", metadata_row_count - nrow(metadata))

message("Discovering Salmon quant files under: ", results_root)
manifest <- discover_quant_files(results_root)

message("Discovering sample-level MultiQC metrics under: ", results_root)
qc_metrics <- discover_sample_qc_metrics(results_root)

expected_samples <- metadata$sample
observed_samples <- manifest$sample

missing_quant <- setdiff(expected_samples, observed_samples)
extra_quant <- setdiff(observed_samples, expected_samples)
ignored_extra_quant <- extra_quant[sample_mouse_ids(extra_quant) %in% excluded_mouse_ids]
extra_quant <- setdiff(extra_quant, ignored_extra_quant)
if (length(ignored_extra_quant) > 0) {
  message("Ignoring Salmon quant files for excluded mice: ", length(ignored_extra_quant))
}
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
ignored_extra_qc <- extra_qc[sample_mouse_ids(extra_qc) %in% excluded_mouse_ids]
extra_qc <- setdiff(extra_qc, ignored_extra_qc)
if (length(ignored_extra_qc) > 0) {
  message("Ignoring MultiQC sample metrics for excluded mice: ", length(ignored_extra_qc))
}
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

# Reorder everything to the metadata sample order once, so all later matrices
# and metadata tables stay aligned.
sample_order <- metadata |>
  dplyr::select("sample")

manifest <- sample_order |>
  dplyr::left_join(manifest, by = "sample")
qc_metrics <- sample_order |>
  dplyr::left_join(qc_metrics, by = "sample")

stopifnot(identical(metadata$sample, manifest$sample))
stopifnot(identical(metadata$sample, qc_metrics$sample))

# Fold technical covariates directly into the cleaned metadata table so later
# DEA scripts only need one metadata object.
metadata <- metadata |>
  dplyr::left_join(qc_metrics, by = "sample")

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
  dropInfReps = TRUE,
  ignoreTxVersion = FALSE
)

# Gene-name annotations are convenient later when writing marker lookups and DE
# result tables.
gene_annotation <- tx2gene |>
  dplyr::distinct(.data$gene_id, .keep_all = TRUE) |>
  dplyr::select("gene_id", "gene_name") |>
  dplyr::arrange(.data$gene_id)

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

population_counts <- metadata |>
  dplyr::count(.data$population, name = "n_samples")
write_tsv(population_counts, file.path(output_dir, "population_sample_counts.tsv"))
writeLines(capture.output(sessionInfo()), con = file.path(output_dir, "sessionInfo.txt"))

message("Finished building unified Salmon gene-level input object")
message("Samples imported: ", nrow(metadata))
message("Genes imported: ", nrow(txi_gene$counts))
message("Batches merged: ", paste(sort(unique(manifest$batch)), collapse = ", "))
