#!/usr/bin/env Rscript

# Figure 02: two-panel summary of OXY ensemble variance.
# Left: contour map of sample variance across ensemble members.
# Right: histogram of grid-cell variances with summary-statistics inlay.

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
  } else {
    normalizePath(getwd())
  }
})

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

source(file.path(script_dir, "plot_utils.R"))

figure_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(file.path(figure_dir, ".."), winslash = "/", mustWork = TRUE)
input_dir <- file.path(project_dir, "Data", "inputs")
ensemble_file <- Sys.getenv(
  "FK_ENSEMBLE_FILE",
  unset = file.path(input_dir, "oxy_ensemble_20_members.rds")
)
if (!file.exists(ensemble_file)) {
  stop("Missing ensemble file: ", ensemble_file, call. = FALSE)
}

ens <- readRDS(ensemble_file)
required_fields <- c("domain", "members", "fields")
if (!all(required_fields %in% names(ens))) {
  stop(
    "Ensemble RDS is missing required fields. Expected: ",
    paste(required_fields, collapse = ", "),
    call. = FALSE
  )
}

if (length(dim(ens$fields)) != 3L) {
  stop("Expected ens$fields to be a 3D array [R x C x M].", call. = FALSE)
}

lat_keep <- as.numeric(ens$domain$lat)
lon_keep <- as.numeric(ens$domain$lon)
R <- length(lat_keep)
C <- length(lon_keep)
M <- length(ens$members)

field_dims <- dim(ens$fields)
if (!identical(field_dims[1:2], c(R, C))) {
  stop("Grid dimensions in ens$fields do not match domain lat/lon lengths.", call. = FALSE)
}
if (field_dims[3] != M) {
  stop("Number of ensemble members does not match third dimension of ens$fields.", call. = FALSE)
}

Xp <- matrix(rep(lon_keep, each = R), nrow = R, ncol = C)
Yp <- matrix(rep(lat_keep, times = C), nrow = R, ncol = C)

plot_lon_limits <- if (!is.null(ens$domain$plot_lon_limits)) {
  as.numeric(ens$domain$plot_lon_limits)
} else {
  range(lon_keep, na.rm = TRUE) + c(-2, 2)
}
plot_lat_limits <- if (!is.null(ens$domain$plot_lat_limits)) {
  as.numeric(ens$domain$plot_lat_limits)
} else {
  range(lat_keep, na.rm = TRUE) + c(-2, 2)
}

world_df <- build_world_boundaries(
  lon_limits = plot_lon_limits,
  lat_limits = plot_lat_limits,
  lon_buffer = 0,
  lat_buffer = 0
)

var_mat <- apply(ens$fields, c(1, 2), stats::var, na.rm = TRUE)
interior_mask <- matrix(TRUE, nrow = R, ncol = C)
interior_mask[c(1, R), ] <- FALSE
interior_mask[, c(1, C)] <- FALSE

var_vals <- as.vector(var_mat[interior_mask])
var_vals <- var_vals[is.finite(var_vals)]

if (length(var_vals) == 0L) {
  stop("No finite variance values were computed.", call. = FALSE)
}

var_limits <- c(0, max(var_vals, na.rm = TRUE))
if (!is.finite(var_limits[2]) || var_limits[2] <= 0) {
  var_limits[2] <- 1
}

variance_palette <- c("#f7fbff", "#d6e6f5", "#92c5de", "#4393c3", "#2166ac", "#053061")

contour_df <- data.frame(
  X = as.vector(Xp),
  Y = as.vector(Yp),
  Z = as.vector(var_mat)
)

contour_plot <- ggplot2::ggplot(contour_df, ggplot2::aes(x = X, y = Y, z = Z)) +
  ggplot2::geom_contour_filled(
    ggplot2::aes(fill = after_stat((level_low + level_high) / 2)),
    bins = 42
  ) +
  ggplot2::coord_equal(xlim = plot_lon_limits, ylim = plot_lat_limits, expand = FALSE) +
  ggplot2::scale_x_continuous(expand = c(0, 0), limits = plot_lon_limits) +
  ggplot2::scale_y_continuous(expand = c(0, 0), limits = plot_lat_limits) +
  ggplot2::scale_fill_gradientn(
    colors = variance_palette,
    limits = var_limits,
    oob = scales::squish
  ) +
  ggplot2::labs(
    title = "(A) Sample variance across ensemble members",
    fill = NULL
  ) +
  ggplot2::theme_bw(base_size = 9, base_family = "serif") +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = "black", linewidth = 0.35),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 8.5, color = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.25),
    plot.title = ggplot2::element_text(hjust = 0),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 8.5),
    legend.text = ggplot2::element_text(size = 7.5),
    legend.key.height = grid::unit(22, "pt"),
    plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "pt")
  )

contour_plot <- add_world_boundaries(contour_plot, world_df = world_df)

q <- stats::quantile(var_vals, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
stats_label <- paste(
  sprintf("mean   = %.3f", mean(var_vals)),
  sprintf("sd     = %.3f", stats::sd(var_vals)),
  sep = "\n"
)

hist_plot <- ggplot2::ggplot(
  data.frame(variance = var_vals),
  ggplot2::aes(x = variance)
) +
  ggplot2::geom_histogram(
    bins = 42,
    fill = "#3f7f93",
    color = "white",
    linewidth = 0.2
  ) +
  ggplot2::annotate(
    geom = "label",
    x = Inf,
    y = Inf,
    hjust = 1.04,
    vjust = 1.04,
    label = stats_label,
    size = 2.95,
    family = "mono",
    label.size = 0.2,
    lineheight = 1.0,
    fill = grDevices::adjustcolor("white", alpha.f = 0.9)
  ) +
  ggplot2::labs(
    title = "(B) Distribution of sample variances",
    x = "Variance (OXY^2)",
    y = "Count"
  ) +
  ggplot2::theme_bw(base_size = 9, base_family = "serif") +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = "black", linewidth = 0.35),
    axis.text = ggplot2::element_text(size = 8.5, color = "black"),
    axis.title = ggplot2::element_text(size = 9),
    axis.ticks = ggplot2::element_line(linewidth = 0.25),
    plot.title = ggplot2::element_text(hjust = 0),
    plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "pt")
  )

figure_width_in <- 10.2
figure_height_in <- 4.8
layout_widths <- c(1.25, 1)

decompose_width <- function(u) {
  type <- grid::unitType(u)
  if (type == "null") {
    return(list(fixed = 0, null = as.numeric(u)))
  }
  if (type == "sum") {
    parts <- unclass(u)[[1]][[2]]
    decomposed <- lapply(parts, decompose_width)
    return(list(
      fixed = sum(vapply(decomposed, `[[`, numeric(1), "fixed")),
      null = sum(vapply(decomposed, `[[`, numeric(1), "null"))
    ))
  }
  list(fixed = grid::convertWidth(u, "in", valueOnly = TRUE), null = 0)
}

compute_panel_width_in <- function(grob, slot_width_in) {
  panel_cols <- unique(grob$layout$l[grob$layout$name == "panel"])
  total_fixed <- 0
  total_null <- 0
  panel_fixed <- 0
  panel_null <- 0

  for (i in seq_along(grob$widths)) {
    dec <- decompose_width(grob$widths[i])
    total_fixed <- total_fixed + dec$fixed
    total_null <- total_null + dec$null
    if (i %in% panel_cols) {
      panel_fixed <- panel_fixed + dec$fixed
      panel_null <- panel_null + dec$null
    }
  }

  if (total_null <= 0) {
    return(panel_fixed)
  }
  panel_fixed + (slot_width_in - total_fixed) * panel_null / total_null
}

slot_width_contour_in <- figure_width_in * layout_widths[1] / sum(layout_widths)
slot_width_hist_in <- figure_width_in * layout_widths[2] / sum(layout_widths)

probe_pdf <- tempfile(fileext = ".pdf")
grDevices::pdf(file = probe_pdf, width = figure_width_in, height = figure_height_in)
contour_probe_grob <- ggplot2::ggplotGrob(contour_plot)
hist_probe_grob <- ggplot2::ggplotGrob(hist_plot)
contour_panel_width_in <- compute_panel_width_in(contour_probe_grob, slot_width_contour_in)
hist_panel_width_in <- compute_panel_width_in(hist_probe_grob, slot_width_hist_in)

map_aspect_ratio <- diff(plot_lat_limits) / diff(plot_lon_limits)
hist_aspect_ratio <- map_aspect_ratio * contour_panel_width_in / hist_panel_width_in

hist_plot <- hist_plot +
  ggplot2::theme(aspect.ratio = hist_aspect_ratio)

figure_grob <- suppressWarnings(gridExtra::arrangeGrob(
  contour_plot,
  hist_plot,
  nrow = 1,
  ncol = 2,
  widths = layout_widths
))
grDevices::dev.off()
unlink(probe_pdf)

output_dir <- Sys.getenv("FK_OUTPUT_DIR", unset = file.path(figure_dir, "out"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_pdf <- file.path(output_dir, "Figure02.pdf")
output_png <- file.path(output_dir, "Figure02.png")

pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
suppressWarnings(ggplot2::ggsave(
  filename = output_pdf,
  plot = figure_grob,
  device = pdf_device,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  bg = "white"
))
crop_pdf_in_place(output_pdf)

suppressWarnings(ggplot2::ggsave(
  filename = output_png,
  plot = figure_grob,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  dpi = 600,
  bg = "white"
))

message("Saved figure PDF: ", normalizePath(output_pdf, mustWork = FALSE))
message("Saved figure PNG: ", normalizePath(output_png, mustWork = FALSE))
