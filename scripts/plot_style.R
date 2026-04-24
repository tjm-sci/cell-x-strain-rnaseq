# Shared plot styling helpers for the EXP383 RNA-seq project.
#
# These settings deliberately mirror the visual conventions used in the
# neighbouring FANS analysis repo:
#   ../exp383-analysis_of_FANS_data
#
# Function:
# - keep population colours stable across the project
# - keep ggplot typography and panel styling consistent
# - make the style easy to reuse from multiple scripts

# Population palette copied from the FANS analysis repo. We keep the names used
# in this RNA-seq project as the primary keys, but also expose a few aliases so
# downstream scripts can stay robust to small naming differences.
exp383_population_palette <- c(
  "NeuN" = "#16a829",
  "PU1" = "#eeaf0f",
  "SOX10" = "#dd0d0d",
  "SOX2" = "#e7298a",
  "NeuN+" = "#16a829",
  "PU1+" = "#eeaf0f",
  "PU.1" = "#eeaf0f",
  "PU.1+" = "#eeaf0f",
  "SOX10+" = "#dd0d0d",
  "SOX2+" = "#e7298a"
)

# Preferred display order for the four sorted nuclei populations used here.
exp383_population_levels <- c("NeuN", "SOX10", "SOX2", "PU1")
exp383_population_palette_clean <- exp383_population_palette[exp383_population_levels]
exp383_population_display_labels <- c(
  "NeuN" = "NeuN+",
  "SOX10" = "SOX10+",
  "SOX2" = "SOX2+",
  "PU1" = "PU1+"
)
exp383_population_display_levels <- unname(exp383_population_display_labels[exp383_population_levels])

# Helper to standardise factor ordering before plotting.
exp383_population_factor <- function(x) {
  factor(as.character(x), levels = exp383_population_levels)
}

exp383_population_display_factor <- function(x) {
  display_values <- unname(exp383_population_display_labels[as.character(x)])
  factor(display_values, levels = exp383_population_display_levels)
}

# Base ggplot theme. This keeps the light-panel look from the FANS repo and the
# same text sizing choices used there, while also removing minor grid clutter.
exp383_theme <- function(base_size = 12) {
  ggplot2::theme_light(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 11, colour = "black"),
      strip.text = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Manual scales for any plot where population is encoded directly.
exp383_scale_colour_population <- function(...) {
  ggplot2::scale_colour_manual(
    values = exp383_population_palette_clean,
    breaks = exp383_population_levels,
    labels = exp383_population_display_labels[exp383_population_levels],
    drop = FALSE,
    ...
  )
}

exp383_scale_fill_population <- function(...) {
  ggplot2::scale_fill_manual(
    values = exp383_population_palette_clean,
    breaks = exp383_population_levels,
    labels = exp383_population_display_labels[exp383_population_levels],
    drop = FALSE,
    ...
  )
}

exp383_scale_x_population <- function(...) {
  ggplot2::scale_x_discrete(
    limits = exp383_population_levels,
    labels = exp383_population_display_labels,
    drop = FALSE,
    ...
  )
}

# pheatmap accepts a named list of annotation colours. We only hard-code the
# population palette here; other annotations can keep their defaults.
exp383_population_annotation_colors <- function() {
  display_palette <- stats::setNames(
    unname(exp383_population_palette_clean),
    exp383_population_display_levels
  )

  list(population = display_palette)
}

# High-resolution raster settings for project figures. PNG is used throughout
# because it is broadly compatible with slides, documents, and web exports.
exp383_plot_dpi <- 300

exp383_save_ggplot <- function(output_file, plot, width, height, dpi = exp383_plot_dpi) {
  ggplot2::ggsave(
    filename = output_file,
    plot = plot,
    device = "png",
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

exp383_open_png_device <- function(output_file, width, height, dpi = exp383_plot_dpi) {
  grDevices::png(
    filename = output_file,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white"
  )
}
