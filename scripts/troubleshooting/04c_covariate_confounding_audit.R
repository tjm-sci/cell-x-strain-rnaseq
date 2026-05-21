#!/usr/bin/env Rscript

# Quantify and visualise associations between retained covariates and biology.

suppressPackageStartupMessages({
  required_packages <- c("dplyr", "ggplot2", "here")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }
})

suppressMessages(here::i_am("scripts/troubleshooting/04c_covariate_confounding_audit.R"))
source(here::here("scripts", "path_helpers.R"))
source(here::here("scripts", "plot_style.R"))
source(here::here("scripts", "03a_deseq2_dea_utilities.R"))
source(here::here("scripts", "troubleshooting", "04a_deseq2_troubleshooting_utilities.R"))

SETTINGS <- list(
  population = "all",
  input_rds = trouble_input_rds_default,
  screening_root = trouble_screening_root_default,
  output_root = file.path(trouble_root_default, "02_covariate_confounding")
)

args <- parse_key_value_args(SETTINGS)
input_rds <- resolve_project_path(args$input_rds, must_work = TRUE)
screening_root <- resolve_project_path(args$screening_root, must_work = TRUE)
output_root <- resolve_project_path(args$output_root)

analysis_input <- load_dea_input(input_rds)
metadata <- add_biological_fields(prepare_coldata_for_deseq(analysis_input$sample_metadata))
available_populations <- order_known_levels(metadata$population, exp383_population_levels)
populations_to_run <- parse_population_selection(args$population, available_populations)

ensure_dir(output_root)

safe_p <- function(expr) {
  tryCatch(expr, error = function(e) NA_real_)
}

numeric_association <- function(data, covariate_value_col, biological_var) {
  form <- stats::as.formula(paste(covariate_value_col, "~", biological_var))
  fit <- stats::lm(form, data = data)
  aov_table <- stats::anova(fit)
  data.frame(
    test_type = "linear_model_anova",
    p_value = safe_p(aov_table$`Pr(>F)`[[1]]),
    effect_size_or_r2 = safe_p(summary(fit)$r.squared),
    notes = "R2 from lm(covariate ~ biological_variable).",
    stringsAsFactors = FALSE
  )
}

categorical_association <- function(data, covariate_col, biological_var) {
  tab <- table(data[[covariate_col]], data[[biological_var]])
  if (any(dim(tab) < 2)) {
    return(data.frame(
      test_type = "not_tested",
      p_value = NA_real_,
      effect_size_or_r2 = NA_real_,
      notes = "Fewer than two levels observed.",
      stringsAsFactors = FALSE
    ))
  }

  chisq <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  fisher_p <- safe_p(stats::fisher.test(tab)$p.value)
  use_fisher <- any(chisq$expected < 5) && !is.na(fisher_p)
  p_value <- if (use_fisher) fisher_p else chisq$p.value
  n <- sum(tab)
  cramer_v <- sqrt(unname(chisq$statistic) / (n * (min(dim(tab)) - 1)))
  data.frame(
    test_type = if (use_fisher) "fisher_exact" else "chi_squared",
    p_value = p_value,
    effect_size_or_r2 = cramer_v,
    notes = "Effect size is Cramer's V.",
    stringsAsFactors = FALSE
  )
}

plot_numeric_covariate <- function(data, covariate_col, covariate_label, biological_var, output_file) {
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[biological_var]], y = .data[[covariate_col]], fill = .data[[biological_var]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.75, show.legend = FALSE) +
    ggplot2::geom_jitter(width = 0.16, height = 0, size = 1.8, alpha = 0.85, show.legend = FALSE) +
    exp383_theme(base_size = 11) +
    ggplot2::scale_fill_manual(values = set3_palette(length(unique(data[[biological_var]])))) +
    ggplot2::labs(
      title = paste(covariate_label, "by", biological_var),
      subtitle = "Retained covariate distribution across biological groups.",
      x = biological_var,
      y = covariate_label
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  exp383_save_ggplot(output_file, plot = p, width = 8, height = 5.5)
}

plot_categorical_covariate <- function(data, covariate_col, covariate_label, biological_var, output_prefix) {
  bar_data <- as.data.frame(table(data[[biological_var]], data[[covariate_col]]), stringsAsFactors = FALSE)
  colnames(bar_data) <- c("biological_value", "covariate_value", "n")
  bar_data <- bar_data[bar_data$n > 0, , drop = FALSE]
  bar_data$biological_value <- factor(bar_data$biological_value, levels = unique(as.character(data[[biological_var]])))

  cov_levels <- unique(as.character(bar_data$covariate_value))
  p_bar <- ggplot2::ggplot(bar_data, ggplot2::aes(x = biological_value, y = n, fill = covariate_value)) +
    ggplot2::geom_col(position = "fill") +
    exp383_theme(base_size = 11) +
    ggplot2::scale_fill_manual(values = stats::setNames(set3_palette(length(cov_levels)), cov_levels)) +
    ggplot2::labs(
      title = paste(covariate_label, "composition by", biological_var),
      subtitle = "Stacked proportions within each biological group.",
      x = biological_var,
      y = "Proportion",
      fill = covariate_label
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  exp383_save_ggplot(paste0(output_prefix, "_stacked_bar.png"), plot = p_bar, width = 8, height = 5.5)

  heat_data <- bar_data
  group_totals <- aggregate(n ~ biological_value, heat_data, sum)
  heat_data$proportion <- heat_data$n / group_totals$n[match(heat_data$biological_value, group_totals$biological_value)]
  p_heat <- ggplot2::ggplot(heat_data, ggplot2::aes(x = biological_value, y = covariate_value, fill = proportion)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = n), size = 3) +
    exp383_theme(base_size = 11) +
    ggplot2::scale_fill_viridis_c(option = "viridis", limits = c(0, 1)) +
    ggplot2::labs(
      title = paste(covariate_label, "contingency heatmap by", biological_var),
      subtitle = "Tile colour is within-group proportion; text is sample count.",
      x = biological_var,
      y = covariate_label,
      fill = "Proportion"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  exp383_save_ggplot(paste0(output_prefix, "_heatmap.png"), plot = p_heat, width = 8, height = 5.5)
}

classify_covariate_risk <- function(covariate_row, association_rows, original_formula) {
  covariate <- covariate_row$covariate
  analysis_column <- covariate_row$analysis_column
  role <- covariate_row$covariate_role
  included <- covariate %in% all.vars(stats::as.formula(original_formula)) ||
    analysis_column %in% all.vars(stats::as.formula(original_formula))

  strong_association <- any(!is.na(association_rows$p_value) & association_rows$p_value < 0.05, na.rm = TRUE)
  max_effect <- suppressWarnings(max(abs(association_rows$effect_size_or_r2), na.rm = TRUE))
  if (!is.finite(max_effect)) max_effect <- NA_real_

  direct_count_like <- grepl("num_mapped|total_sequences|passed_filter_reads", covariate)
  qc_like <- role == "technical_qc" || covariate %in% trouble_qc_covariates
  known_batch <- covariate %in% c("inoculation_batch", "date_nuc_prep_days")

  if (direct_count_like) {
    risk <- "do_not_include_without_strong_justification"
    reason <- "Direct library-size or mapped-read covariate can absorb real abundance shifts."
    handling <- "Avoid in primary DEA unless a documented sensitivity analysis supports inclusion."
  } else if (qc_like && strong_association) {
    risk <- "high_risk"
    reason <- "QC covariate is associated with biological design."
    handling <- "Use only as sensitivity analysis or after explicit sample-quality review."
  } else if (qc_like) {
    risk <- "high_risk"
    reason <- "Sequencing/QC covariate may absorb disease-linked technical or biological signal."
    handling <- "Prefer exclusion/QC sensitivity over routine adjustment."
  } else if (known_batch && strong_association) {
    risk <- "moderate_risk"
    reason <- "Known experimental variable is at least partly associated with biology."
    handling <- "Consider in safer selected design, but compare against minimal model."
  } else if (known_batch) {
    risk <- "low_risk"
    reason <- "Known experimental variable without strong biological association."
    handling <- "Reasonable to include if it improves residual structure."
  } else if (strong_association || (!is.na(max_effect) && max_effect >= 0.2)) {
    risk <- "high_risk"
    reason <- "Covariate shows association with biological design."
    handling <- "Avoid without strong biological and technical justification."
  } else {
    risk <- "low_risk"
    reason <- "No strong biological association detected in this audit."
    handling <- "Acceptable candidate for reduced sensitivity designs."
  }

  data.frame(
    covariate = covariate,
    retained_in_original_screen = covariate_row$selected,
    included_in_original_formula = included,
    risk_class = risk,
    reason = reason,
    recommended_handling = handling,
    stringsAsFactors = FALSE
  )
}

all_associations <- list()
all_risks <- list()

for (population in populations_to_run) {
  message("Auditing covariate confounding for ", population)
  population_root <- ensure_dir(file.path(output_root, population))
  table_dir <- ensure_dir(file.path(population_root, "tables"))
  plot_dir <- ensure_dir(file.path(population_root, "plots"))

  selected <- read_selected_covariates(screening_root, population)
  selected <- selected[selected$selected %in% c(TRUE, "TRUE", "true", "1", 1), , drop = FALSE]
  original_formula <- read_population_design_formula(screening_root, population)
  coldata <- droplevels(metadata[metadata$population == population, , drop = FALSE])
  coldata <- ensure_formula_z_columns(coldata, original_formula)

  population_associations <- list()
  for (i in seq_len(nrow(selected))) {
    covariate <- selected$covariate[[i]]
    analysis_column <- selected$analysis_column[[i]]
    cov_type <- selected$covariate_type[[i]]
    value_col <- if (covariate %in% colnames(coldata)) covariate else analysis_column
    if (!value_col %in% colnames(coldata)) next

    for (bio_var in intersect(trouble_biological_vars, colnames(coldata))) {
      test_data <- coldata[, c(value_col, bio_var), drop = FALSE]
      test_data <- test_data[stats::complete.cases(test_data), , drop = FALSE]
      if (nrow(test_data) == 0 || length(unique(test_data[[bio_var]])) < 2) next

      if (identical(cov_type, "numeric")) {
        test_data[[value_col]] <- as.numeric(test_data[[value_col]])
        test_result <- numeric_association(test_data, value_col, bio_var)
        plot_file <- file.path(plot_dir, paste0(sanitize_filename(covariate), "_by_", sanitize_filename(bio_var), ".png"))
        plot_numeric_covariate(test_data, value_col, covariate, bio_var, plot_file)
      } else {
        test_result <- categorical_association(test_data, value_col, bio_var)
        output_prefix <- file.path(plot_dir, paste0(sanitize_filename(covariate), "_by_", sanitize_filename(bio_var)))
        plot_categorical_covariate(test_data, value_col, covariate, bio_var, output_prefix)
      }

      population_associations[[paste(covariate, bio_var, sep = "__")]] <- data.frame(
        population = population,
        covariate = covariate,
        analysis_column = analysis_column,
        covariate_type = cov_type,
        biological_variable = bio_var,
        test_type = test_result$test_type,
        p_value = test_result$p_value,
        effect_size_or_r2 = test_result$effect_size_or_r2,
        notes = test_result$notes,
        stringsAsFactors = FALSE
      )
    }
  }

  population_association_table <- if (length(population_associations) > 0) {
    do.call(rbind, population_associations)
  } else {
    data.frame()
  }
  if (nrow(population_association_table) > 0) {
    rownames(population_association_table) <- NULL
    write_tsv(population_association_table, file.path(table_dir, "covariate_biology_association_summary.tsv"))
    all_associations[[population]] <- population_association_table
  }

  population_risks <- list()
  for (i in seq_len(nrow(selected))) {
    covariate <- selected$covariate[[i]]
    rows <- population_association_table[population_association_table$covariate == covariate, , drop = FALSE]
    risk <- classify_covariate_risk(selected[i, , drop = FALSE], rows, original_formula)
    risk$population <- population
    risk$analysis_column <- selected$analysis_column[[i]]
    risk <- risk[, c("population", "covariate", "analysis_column", "retained_in_original_screen", "included_in_original_formula", "risk_class", "reason", "recommended_handling")]
    population_risks[[covariate]] <- risk
  }
  population_risk_table <- do.call(rbind, population_risks)
  rownames(population_risk_table) <- NULL
  write_tsv(population_risk_table, file.path(table_dir, "covariate_risk_classification.tsv"))
  all_risks[[population]] <- population_risk_table
}

if (length(all_associations) > 0) {
  association_summary <- do.call(rbind, all_associations)
  rownames(association_summary) <- NULL
  write_tsv(association_summary, file.path(output_root, "covariate_biology_association_summary.tsv"))
}

if (length(all_risks) > 0) {
  risk_summary <- do.call(rbind, all_risks)
  rownames(risk_summary) <- NULL
  write_tsv(risk_summary, file.path(output_root, "covariate_risk_classification.tsv"))
}

message("Finished covariate confounding audit.")
