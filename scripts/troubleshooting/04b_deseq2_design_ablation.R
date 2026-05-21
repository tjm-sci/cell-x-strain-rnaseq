#!/usr/bin/env Rscript

# Run a DESeq2 design ladder for each population.
#
# This is the first troubleshooting pass: it tests whether signal is present in
# simpler designs and lost after the covariate-screened design is applied.

suppressPackageStartupMessages({
  required_packages <- c("dplyr", "DESeq2", "SummarizedExperiment", "ggplot2", "ggrepel", "pheatmap", "here")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }
})

suppressMessages(here::i_am("scripts/troubleshooting/04b_deseq2_design_ablation.R"))
source(here::here("scripts", "path_helpers.R"))
source(here::here("scripts", "plot_style.R"))
source(here::here("scripts", "03a_deseq2_dea_utilities.R"))
source(here::here("scripts", "troubleshooting", "04a_deseq2_troubleshooting_utilities.R"))

SETTINGS <- list(
  population = "all",
  input_rds = trouble_input_rds_default,
  screening_root = trouble_screening_root_default,
  output_root = file.path(trouble_root_default, "01_deseq2_design_ablation"),
  contrast_tsv = trouble_contrast_tsv_default,
  min_count = as.character(trouble_min_count_default),
  min_samples = "auto",
  padj_cutoff = as.character(trouble_padj_cutoff_default),
  lfc_threshold = as.character(trouble_lfc_threshold_default),
  label_top_n = "20"
)

args <- parse_key_value_args(SETTINGS)

input_rds <- resolve_project_path(args$input_rds, must_work = TRUE)
screening_root <- resolve_project_path(args$screening_root, must_work = TRUE)
output_root <- resolve_project_path(args$output_root)
contrast_tsv <- resolve_project_path(args$contrast_tsv, must_work = TRUE)

min_count <- as.integer(args$min_count)
min_samples <- if (identical(args$min_samples, "auto")) "auto" else as.integer(args$min_samples)
padj_cutoff <- as.numeric(args$padj_cutoff)
lfc_threshold <- as.numeric(args$lfc_threshold)
label_top_n <- as.integer(args$label_top_n)

analysis_input <- load_dea_input(input_rds)
metadata <- add_biological_fields(prepare_coldata_for_deseq(analysis_input$sample_metadata))
available_populations <- order_known_levels(metadata$population, exp383_population_levels)
populations_to_run <- parse_population_selection(args$population, available_populations)
contrast_table <- read_contrast_table(contrast_tsv)

ensure_dir(output_root)

all_model_summaries <- list()
all_contrast_summaries <- list()
all_design_tables <- list()
skipped_models <- list()

for (population in populations_to_run) {
  message("Running design ablation for ", population)
  coldata_pop <- droplevels(metadata[metadata$population == population, , drop = FALSE])
  design_table <- build_design_ladder(screening_root, population, coldata_pop)
  all_design_tables[[population]] <- cbind(population = population, design_table)

  population_root <- ensure_dir(file.path(output_root, population))
  write_tsv(design_table, file.path(population_root, "design_ladder.tsv"))

  for (i in seq_len(nrow(design_table))) {
    design_id <- design_table$design_id[[i]]
    design_label <- design_table$design_label[[i]]
    model_formula <- design_table$model_formula[[i]]
    design_root <- ensure_dir(file.path(population_root, design_id))

    message("  ", population, " / ", design_id)
    result <- tryCatch(
      {
        prepared <- prepare_population_troubleshooting(
          analysis_input = analysis_input,
          population = population,
          design_formula = model_formula,
          min_count = min_count,
          min_samples = min_samples
        )
        model_summary <- fit_troubleshooting_model(
          prepared = prepared,
          design_formula = model_formula,
          population = population,
          output_dir = design_root,
          design_id = design_id,
          design_label = design_label
        )
        dds <- readRDS(file.path(design_root, "model_fit", "dds.rds"))
        contrast_summary <- export_troubleshooting_contrasts(
          dds = dds,
          analysis_input = analysis_input,
          contrast_table = contrast_table,
          output_dir = design_root,
          population = population,
          design_id = design_id,
          design_label = design_label,
          model_formula = model_formula,
          padj_cutoff = padj_cutoff,
          lfc_threshold = lfc_threshold,
          label_top_n = label_top_n
        )
        list(model_summary = model_summary, contrast_summary = contrast_summary, error = NULL)
      },
      error = function(e) {
        list(model_summary = NULL, contrast_summary = NULL, error = conditionMessage(e))
      }
    )

    if (!is.null(result$error)) {
      warning(sprintf("Skipping %s %s: %s", population, design_id, result$error), call. = FALSE)
      skipped_models[[paste(population, design_id, sep = "__")]] <- data.frame(
        population = population,
        design_id = design_id,
        design_label = design_label,
        model_formula = model_formula,
        skip_reason = result$error,
        stringsAsFactors = FALSE
      )
      next
    }

    all_model_summaries[[paste(population, design_id, sep = "__")]] <- result$model_summary
    all_contrast_summaries[[paste(population, design_id, sep = "__")]] <- result$contrast_summary
  }
}

if (length(all_design_tables) > 0) {
  design_table_all <- do.call(rbind, all_design_tables)
  rownames(design_table_all) <- NULL
  write_tsv(design_table_all, file.path(output_root, "design_ladder_all_populations.tsv"))
}

if (length(all_model_summaries) > 0) {
  model_summary_all <- do.call(rbind, all_model_summaries)
  rownames(model_summary_all) <- NULL
  write_tsv(model_summary_all, file.path(output_root, "model_fit_summary.tsv"))
}

if (length(all_contrast_summaries) > 0) {
  contrast_summary_all <- do.call(rbind, all_contrast_summaries)
  rownames(contrast_summary_all) <- NULL
  write_tsv(contrast_summary_all, file.path(output_root, "design_ablation_contrast_summary.tsv"))
}

if (length(skipped_models) > 0) {
  skipped <- do.call(rbind, skipped_models)
  rownames(skipped) <- NULL
  write_tsv(skipped, file.path(output_root, "skipped_models.tsv"))
}

message("Finished design ablation.")
