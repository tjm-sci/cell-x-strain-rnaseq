#!/usr/bin/env Rscript

# Write DESeq2 result-state diagnostics, the interim troubleshooting summary,
# and STAR-Salmon planning notes.

suppressPackageStartupMessages({
  required_packages <- c("dplyr", "DESeq2", "SummarizedExperiment", "here")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }
})

suppressMessages(here::i_am("scripts/troubleshooting/04g_write_troubleshooting_summary.R"))
source(here::here("scripts", "path_helpers.R"))
source(here::here("scripts", "plot_style.R"))
source(here::here("scripts", "03a_deseq2_dea_utilities.R"))
source(here::here("scripts", "troubleshooting", "04a_deseq2_troubleshooting_utilities.R"))

SETTINGS <- list(
  troubleshooting_root = trouble_root_default,
  design_ablation_root = file.path(trouble_root_default, "01_deseq2_design_ablation"),
  confounding_root = file.path(trouble_root_default, "02_covariate_confounding"),
  sample_qc_root = file.path(trouble_root_default, "03_sample_qc_sensitivity"),
  combined_root = file.path(trouble_root_default, "04_combined_population_model"),
  positive_control_root = file.path(trouble_root_default, "05_positive_control_gene_audit"),
  output_root = file.path(trouble_root_default, "06_summary_reports"),
  star_plan_root = file.path(trouble_root_default, "99_nfcore_star_salmon_planning"),
  contrast_tsv = trouble_contrast_tsv_default,
  padj_cutoff = as.character(trouble_padj_cutoff_default)
)

args <- parse_key_value_args(SETTINGS)
troubleshooting_root <- resolve_project_path(args$troubleshooting_root)
design_ablation_root <- resolve_project_path(args$design_ablation_root, must_work = TRUE)
confounding_root <- resolve_project_path(args$confounding_root)
sample_qc_root <- resolve_project_path(args$sample_qc_root)
combined_root <- resolve_project_path(args$combined_root)
positive_control_root <- resolve_project_path(args$positive_control_root)
output_root <- resolve_project_path(args$output_root)
star_plan_root <- resolve_project_path(args$star_plan_root)
contrast_tsv <- resolve_project_path(args$contrast_tsv, must_work = TRUE)
padj_cutoff <- as.numeric(args$padj_cutoff)

ensure_dir(output_root)
ensure_dir(star_plan_root)
contrast_table <- read_contrast_table(contrast_tsv)

diagnostic_rows <- list()

for (population in exp383_population_levels) {
  design_ladder_file <- file.path(design_ablation_root, population, "design_ladder.tsv")
  if (!file.exists(design_ladder_file)) next
  design_ladder <- utils::read.delim(design_ladder_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

  for (i in seq_len(nrow(design_ladder))) {
    design_id <- design_ladder$design_id[[i]]
    dds_file <- file.path(design_ablation_root, population, design_id, "model_fit", "dds.rds")
    if (!file.exists(dds_file)) next
    dds <- readRDS(dds_file)
    group_levels <- levels(droplevels(SummarizedExperiment::colData(dds)$group_assignment))

    for (j in seq_len(nrow(contrast_table))) {
      contrast_row <- contrast_table[j, , drop = FALSE]
      if (!all(c(contrast_row$numerator_group, contrast_row$denominator_group) %in% group_levels)) next

      res <- DESeq2::results(
        dds,
        contrast = c("group_assignment", contrast_row$numerator_group, contrast_row$denominator_group),
        alpha = padj_cutoff
      )
      diagnostic_res <- DESeq2::results(
        dds,
        contrast = c("group_assignment", contrast_row$numerator_group, contrast_row$denominator_group),
        independentFiltering = FALSE,
        cooksCutoff = FALSE,
        alpha = padj_cutoff
      )

      diagnostic_rows[[paste(population, design_id, contrast_row$contrast_id, sep = "__")]] <- data.frame(
        population = population,
        design_id = design_id,
        contrast = contrast_row$contrast_id,
        n_total_rows = length(res$pvalue),
        n_all_zero_or_not_tested = sum(is.na(res$pvalue)),
        n_pvalue_na = sum(is.na(res$pvalue)),
        n_padj_na = sum(is.na(res$padj)),
        n_padj_na_but_pvalue_present = sum(!is.na(res$pvalue) & is.na(res$padj)),
        n_nominal_p_lt_0_05 = sum(!is.na(res$pvalue) & res$pvalue < 0.05),
        n_padj_lt_0_05 = sum(!is.na(res$padj) & res$padj < 0.05),
        n_nominal_p_lt_0_05_diagnostic_no_filter_no_cooks = sum(!is.na(diagnostic_res$pvalue) & diagnostic_res$pvalue < 0.05),
        n_padj_lt_0_05_diagnostic_no_filter_no_cooks = sum(!is.na(diagnostic_res$padj) & diagnostic_res$padj < 0.05),
        notes = "Diagnostic result disables independent filtering and Cook's cutoff.",
        stringsAsFactors = FALSE
      )
    }
  }
}

result_state <- do.call(rbind, diagnostic_rows)
rownames(result_state) <- NULL
write_tsv(result_state, file.path(output_root, "deseq2_result_state_audit.tsv"))

read_optional <- function(path) {
  if (!file.exists(path)) return(data.frame())
  utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

design_summary <- read_optional(file.path(design_ablation_root, "design_ablation_contrast_summary.tsv"))
risk_summary <- read_optional(file.path(confounding_root, "covariate_risk_classification.tsv"))
sample_summary <- read_optional(file.path(sample_qc_root, "sample_qc_sensitivity_contrast_summary.tsv"))
combined_summary <- read_optional(file.path(combined_root, "combined_population_contrast_summary.tsv"))
positive_summary <- read_optional(file.path(positive_control_root, "positive_control_summary_by_design.tsv"))

classify_design_ablation_patterns <- function(design_summary) {
  if (nrow(design_summary) == 0) return(data.frame())

  split_rows <- split(design_summary, paste(design_summary$population, design_summary$contrast, sep = "__"))
  pattern_rows <- lapply(split_rows, function(x) {
    population <- x$population[[1]]
    contrast <- x$contrast[[1]]
    model_hits <- stats::setNames(x$n_padj_lt_0_05, x$design_id)
    model_nominal <- stats::setNames(x$n_nominal_p_lt_0_05, x$design_id)

    get_value <- function(named_vector, name) {
      if (name %in% names(named_vector)) return(named_vector[[name]])
      0
    }

    minimal_hits <- get_value(model_hits, "design_01_minimal_group")
    screened_hits <- get_value(model_hits, "design_06_original_screened_design")
    max_hits <- max(model_hits, na.rm = TRUE)
    max_nominal <- max(model_nominal, na.rm = TRUE)
    batch_date_non_qc_hits <- max(model_hits[names(model_hits) %in% c(
      "design_02_known_batch_group",
      "design_03_date_group",
      "design_04_selected_non_qc_covariates_group"
    )], na.rm = TRUE)
    if (!is.finite(batch_date_non_qc_hits)) batch_date_non_qc_hits <- 0

    if (minimal_hits > 0 && screened_hits < minimal_hits) {
      pattern_code <- "A"
      pattern_label <- "Signal present in minimal model but reduced in screened model"
    } else if (max_hits == 0 && max_nominal == 0) {
      pattern_code <- "B"
      pattern_label <- "Signal absent in all tested models"
    } else if (minimal_hits == 0 && batch_date_non_qc_hits > 0) {
      pattern_code <- "C"
      pattern_label <- "Signal only appears after batch/date/non-QC adjustment"
    } else if (max_hits == 0 && max_nominal > 0) {
      pattern_code <- "D"
      pattern_label <- "Nominal signal present but weak after FDR correction"
    } else {
      pattern_code <- "E"
      pattern_label <- "FDR signal is contrast-specific or retained across adjusted models"
    }

    data.frame(
      population = population,
      contrast = contrast,
      minimal_n_padj_lt_0_05 = minimal_hits,
      screened_n_padj_lt_0_05 = screened_hits,
      max_model_n_padj_lt_0_05 = max_hits,
      max_model_n_nominal_p_lt_0_05 = max_nominal,
      batch_date_non_qc_max_n_padj_lt_0_05 = batch_date_non_qc_hits,
      pattern_code = pattern_code,
      pattern_label = pattern_label,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, pattern_rows)
  rownames(out) <- NULL
  out[order(out$population, out$contrast), , drop = FALSE]
}

design_patterns <- classify_design_ablation_patterns(design_summary)
write_tsv(design_patterns, file.path(output_root, "design_ablation_signal_patterns.tsv"))

sum_by_design <- if (nrow(design_summary) > 0) {
  aggregate(n_padj_lt_0_05 ~ population + design_id, design_summary, sum)
} else {
  data.frame()
}

minimal_total <- if (nrow(sum_by_design) > 0) sum(sum_by_design$n_padj_lt_0_05[sum_by_design$design_id == "design_01_minimal_group"]) else NA_integer_
screened_total <- if (nrow(sum_by_design) > 0) sum(sum_by_design$n_padj_lt_0_05[sum_by_design$design_id == "design_06_original_screened_design"]) else NA_integer_

high_risk_covariates <- if (nrow(risk_summary) > 0) {
  risk_summary[risk_summary$risk_class %in% c("high_risk", "do_not_include_without_strong_justification"), , drop = FALSE]
} else {
  data.frame()
}

combined_delta <- if (nrow(combined_summary) > 0) {
  sum(combined_summary$combined_model_n_padj_lt_0_05, na.rm = TRUE) -
    sum(combined_summary$split_original_n_padj_lt_0_05, na.rm = TRUE)
} else {
  NA_real_
}

positive_fdr <- if (nrow(positive_summary) > 0) sum(positive_summary$n_fdr_significant, na.rm = TRUE) else NA_integer_
positive_nominal <- if (nrow(positive_summary) > 0) sum(positive_summary$n_nominal_only, na.rm = TRUE) else NA_integer_
state_filter_delta <- if (nrow(result_state) > 0) {
  sum(result_state$n_padj_lt_0_05_diagnostic_no_filter_no_cooks, na.rm = TRUE) -
    sum(result_state$n_padj_lt_0_05, na.rm = TRUE)
} else {
  NA_real_
}

summary_lines <- c(
  "# EXP383 DESeq2 troubleshooting interim summary",
  "",
  "## Scope",
  "",
  "This summary covers the DESeq2-first troubleshooting pass under `results/dea/troubleshooting/`. Primary DEA outputs under `results/dea/03_deseq2_dea/` were not used as write targets.",
  "",
  "## Design ablation",
  "",
  sprintf("- Minimal-design total zero-effect FDR hits across all tested contrasts: `%s`.", minimal_total),
  sprintf("- Original screened-design total zero-effect FDR hits across all tested contrasts: `%s`.", screened_total),
  sprintf("- Contrast-level signal-pattern classifications written: `%s`.", file.path(output_root, "design_ablation_signal_patterns.tsv")),
  if (!is.na(minimal_total) && !is.na(screened_total) && minimal_total > screened_total) {
    "- Pattern: minimal models recover more signal than screened models, supporting possible over-adjustment."
  } else if (!is.na(minimal_total) && !is.na(screened_total) && minimal_total < screened_total) {
    "- Pattern: screened models recover more signal than minimal models, arguing against simple global over-adjustment."
  } else {
    "- Pattern: minimal and screened models recover similar total signal."
  },
  "",
  "## Covariate confounding",
  "",
  sprintf("- High-risk or do-not-include covariate instances: `%s`.", nrow(high_risk_covariates)),
  "- Direct count/depth covariates and QC covariates are treated as high-risk in the classification table because they can absorb abundance shifts or disease-linked quality changes.",
  "",
  "## Sample QC sensitivity",
  "",
  if (nrow(sample_summary) > 0) {
    sprintf("- QC sensitivity fits completed for `%s` contrast rows.", nrow(sample_summary))
  } else {
    "- QC sensitivity summary was not available when this report was written."
  },
  "- Interpret sample exclusion only through the objective thresholds in `03_sample_qc_sensitivity/exclusion_sets/exclusion_thresholds.tsv`.",
  "",
  "## Combined-population model",
  "",
  sprintf("- Combined model FDR-hit delta versus split original totals: `%s`.", combined_delta),
  "- The combined model is a sensitivity analysis only; it does not replace population-specific interpretation without review.",
  "",
  "## Positive-control genes",
  "",
  sprintf("- Positive-control FDR-significant audit rows: `%s`.", positive_fdr),
  sprintf("- Positive-control nominal-only audit rows: `%s`.", positive_nominal),
  "- Review `05_positive_control_gene_audit/PU1_positive_control_gene_audit.tsv` to distinguish absent, filtered, low-count, high-SE, and directionally inconsistent genes.",
  "",
  "## DESeq2 result-state audit",
  "",
  sprintf("- Extra FDR hits gained by disabling independent filtering and Cook's cutoff across all diagnostic rows: `%s`.", state_filter_delta),
  "- This diagnostic should not be used as final inference; it only shows whether result-state rules hide genes.",
  "",
  "## Current interpretation",
  "",
  "If design ablation, QC sensitivity, and combined modelling fail to recover broad expected signal, the evidence points away from a simple DESeq2 scripting bug and toward upstream quantification or nuclear-RNA biology. Given transcriptome-only Salmon and low Salmon mapping rates in FANS nuclei, STAR/genome-aligned exon-plus-intron counting remains the main follow-up."
)

writeLines(summary_lines, file.path(output_root, "deseq2_troubleshooting_interim_summary.md"))

storage_audit <- data.frame(
  population = "PU1",
  n_samples = NA_integer_,
  work_dir_peak_gb = NA_real_,
  final_outdir_gb = NA_real_,
  save_align_intermeds = FALSE,
  runtime = NA_character_,
  mean_star_uniquely_mapped_percent = NA_real_,
  mean_salmon_mapping_percent = NA_real_,
  notes = "Planning placeholder; populate after the STAR-Salmon pilot is run.",
  stringsAsFactors = FALSE
)
write_tsv(storage_audit, file.path(star_plan_root, "star_salmon_storage_audit.tsv"))

star_plan <- c(
  "# EXP383 STAR-Salmon troubleshooting rerun plan",
  "",
  "## Purpose",
  "",
  "The existing nf-core run used transcriptome-only Salmon quantification. Because these are FANS-sorted nuclei, mature-transcript-only quantification may miss intronic/unspliced nuclear RNA. The next troubleshooting stage should test a genome-aligned STAR-Salmon or STAR plus featureCounts strategy without overwriting current nf-core outputs.",
  "",
  "## Existing inputs",
  "",
  "- `scripts/00_run_nfcore_rnaseq.sh`",
  "- `params/nfcore_rnaseq_fastp_salmon_rrna.yaml`",
  "- `metadata/batches/exp383_nfcore_rnaseq_PU1.csv`",
  "- `metadata/batches/exp383_nfcore_rnaseq_NeuN.csv`",
  "- `metadata/batches/exp383_nfcore_rnaseq_SOX10.csv`",
  "- `metadata/batches/exp383_nfcore_rnaseq_SOX2.csv`",
  "- `resources/reference/Mus_musculus.GRCm39.115.gtf.gz`",
  "- `resources/reference/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz`",
  "",
  "## Proposed new files and directories",
  "",
  "- Parameters: `params/nfcore_rnaseq_star_salmon_troubleshooting.yaml`",
  "- Logs: `logs/nextflow/troubleshooting_star_salmon/`",
  "- Results: `results/nfcore_rnaseq_star_salmon_troubleshooting/`",
  "- Work directory: `/media/tmurphy/2TB_SSD/exp383/nfcore_rnaseq_star_salmon_work/`",
  "",
  "## Storage-aware execution order",
  "",
  "1. Run a tiny PU1 smoke test with 4-8 samples.",
  "2. If storage is acceptable, run all PU1 samples.",
  "3. Review work-dir peak size, final outdir size, STAR mapping, Salmon mapping, and positive-control gene recovery.",
  "4. Only then consider NeuN, SOX10, and SOX2 one population at a time.",
  "",
  "## Pilot command sketch",
  "",
  "```bash",
  "NXF_WORK=/media/tmurphy/2TB_SSD/exp383/nfcore_rnaseq_star_salmon_work/PU1 \\",
  "mamba run -n rnaseq nextflow \\",
  "  -log logs/nextflow/troubleshooting_star_salmon/PU1.nextflow.log \\",
  "  run nf-core/rnaseq \\",
  "  -r 3.23.0 \\",
  "  -profile conda \\",
  "  -c nextflow.config \\",
  "  -params-file params/nfcore_rnaseq_star_salmon_troubleshooting.yaml \\",
  "  --input metadata/batches/exp383_nfcore_rnaseq_PU1.csv \\",
  "  --outdir results/nfcore_rnaseq_star_salmon_troubleshooting/PU1 \\",
  "  -ansi-log false",
  "```",
  "",
  "## Counting strategy to decide before implementation",
  "",
  "- STAR-Salmon transcript-derived gene counts: useful for comparability, but still transcript-annotation dependent.",
  "- featureCounts exon-only counts: standard gene-level comparator.",
  "- featureCounts gene-body or exon+intron counts: likely critical for nuclear RNA.",
  "- Separate exon-only and exon+intron matrices should be retained so the quantification source is explicit.",
  "",
  "## tximport requirements",
  "",
  "- If Salmon transcript-level output is imported, record `txi$countsFromAbundance`.",
  "- If original Salmon estimated counts are used, construct DESeq2 objects with `DESeqDataSetFromTximport()` so length offsets are retained.",
  "- If `countsFromAbundance = \"lengthScaledTPM\"` or another scaled mode is used, document why and do not mix that mode with primary outputs without clear directory names.",
  "",
  "## Post-rerun comparison",
  "",
  "Compare existing transcriptome-only Salmon/tximport counts against STAR/featureCounts exon-only, STAR/featureCounts exon+intron or gene-body counts, and STAR-Salmon transcript-derived gene counts if retained. For each source, rerun minimal DESeq2, safer selected DESeq2, and the positive-control audit."
)

writeLines(star_plan, file.path(star_plan_root, "nfcore_star_salmon_rerun_plan.md"))

message("Finished troubleshooting summary and planning outputs.")
