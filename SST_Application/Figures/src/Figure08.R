#!/usr/bin/env Rscript

# Figure 08: two-panel summary of ensemble variance.
# Left: contour map of sample variance across ensemble members.
# Right: histogram of grid-cell variances with summary-statistics inlay.

resolve_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) return(NULL)
  normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)
}

resolve_project_root <- function() {
  script_path <- resolve_script_path()
  candidates <- c(
    if (!is.null(script_path)) normalizePath(file.path(dirname(script_path), "..", "..", ".."), mustWork = FALSE) else NA_character_,
    if (!is.null(script_path)) normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE) else NA_character_,
    normalizePath(getwd(), mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
  )
  candidates <- unique(candidates[!is.na(candidates)])

  for (cand in candidates) {
    if (file.exists(file.path(cand, "modules", "plot_sst_maps.R"))) {
      return(cand)
    }
  }
  stop("Could not locate project root containing modules/plot_sst_maps.R")
}

project_root <- resolve_project_root()
setwd(project_root)

project_lib <- file.path(project_root, ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

source(file.path("modules", "plot_utils.R"))
source(file.path("modules", "plot_sst_maps.R"))
source(file.path(project_root, "Figures", "src", "figure_io_utils.R"))

figure_root <- file.path(project_root, "Figures")
fig_dir <- resolve_figure_output_dir(figure_root, run_label = "static")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ensemble_rds <- file.path("Data", "processed", "pacific_sst_ensemble_subindexed_201201.rds")
if (!file.exists(ensemble_rds)) {
  stop("Missing ensemble data file: ", ensemble_rds)
}

ens <- readRDS(ensemble_rds)
required_fields <- c("domain", "members", "fields")
if (!all(required_fields %in% names(ens))) {
  stop(
    "Ensemble RDS is missing required fields. Expected: ",
    paste(required_fields, collapse = ", ")
  )
}

if (length(dim(ens$fields)) != 3L) {
  stop("Expected ens$fields to be a 3D array [R x C x M].")
}

lat_keep <- as.numeric(ens$domain$lat)
lon_keep <- as.numeric(ens$domain$lon)
R <- length(lat_keep)
C <- length(lon_keep)
M <- length(ens$members)

field_dims <- dim(ens$fields)
if (!identical(field_dims[1:2], c(R, C))) {
  stop("Grid dimensions in ens$fields do not match domain lat/lon lengths.")
}
if (field_dims[3] != M) {
  stop("Number of ensemble members does not match third dimension of ens$fields.")
}

Xp <- matrix(rep(lon_keep, each = R), nrow = R, ncol = C)
Yp <- matrix(rep(lat_keep, times = C), nrow = R, ncol = C)

land_file <- file.path("Data", "processed", "pacific_land_mask.csv")
if (file.exists(land_file)) {
  land_raw <- utils::read.csv(land_file, stringsAsFactors = FALSE)
  required_land_cols <- c("s", "lat", "lon", "is_land")
  if (!all(required_land_cols %in% names(land_raw))) {
    stop(
      "Land-mask file is missing required fields. Expected: ",
      paste(required_land_cols, collapse = ", ")
    )
  }
  land_raw <- land_raw[order(land_raw$s), required_land_cols]
  if (!identical(as.integer(land_raw$s), seq_len(R * C))) {
    stop("Land-mask site ids do not align with the ensemble grid indexing.")
  }
  land_df <- data.frame(
    X = as.numeric(land_raw$lon),
    Y = as.numeric(land_raw$lat),
    is_land = as.logical(land_raw$is_land)
  )
} else {
  land_df <- build_land_mask(Xp = Xp, Yp = Yp)
}
ocean_mask <- !as.logical(land_df$is_land)
if (length(ocean_mask) != R * C) {
  stop("Land mask length does not match the ensemble grid size.")
}

world_df <- build_world_boundaries(
  lon_limits = range(Xp, na.rm = TRUE),
  lat_limits = range(Yp, na.rm = TRUE)
)

var_mat <- apply(ens$fields, c(1, 2), stats::var, na.rm = TRUE)
var_limits <- c(0, 4)
var_vals <- as.vector(var_mat)
var_vals <- var_vals[ocean_mask & is.finite(var_vals)]

if (length(var_vals) == 0L) {
  stop("No finite variance values were computed.")
}

variance_palette <- c("#f7fbff", "#d6e6f5", "#92c5de", "#4393c3", "#2166ac", "#053061")

contour_plot <- suppressWarnings(plot_sst_contour_map(
  Z = var_mat,
  Xp = Xp,
  Yp = Yp,
  bins = 42,
  palette = variance_palette,
  limits = var_limits,
  title = "(A) Sample variance across ensemble members",
  land_df = land_df,
  world_df = world_df
)) +
  ggplot2::labs(fill = NULL) +
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
  ggplot2::coord_cartesian(xlim = c(0, 5)) +
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
    x = NULL,
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

contour_probe_grob <- ggplot2::ggplotGrob(contour_plot)
hist_probe_grob <- ggplot2::ggplotGrob(hist_plot)

slot_width_contour_in <- figure_width_in * layout_widths[1] / sum(layout_widths)
slot_width_hist_in <- figure_width_in * layout_widths[2] / sum(layout_widths)

grDevices::pdf(NULL, width = figure_width_in, height = figure_height_in)
contour_panel_width_in <- compute_panel_width_in(contour_probe_grob, slot_width_contour_in)
hist_panel_width_in <- compute_panel_width_in(hist_probe_grob, slot_width_hist_in)
grDevices::dev.off()

map_aspect_ratio <- diff(range(Yp, na.rm = TRUE)) / diff(range(Xp, na.rm = TRUE))
hist_aspect_ratio <- 1.08 * map_aspect_ratio * contour_panel_width_in / hist_panel_width_in

hist_plot <- hist_plot +
  ggplot2::theme(aspect.ratio = hist_aspect_ratio)

figure_grob <- suppressWarnings(gridExtra::arrangeGrob(
  contour_plot,
  hist_plot,
  nrow = 1,
  ncol = 2,
  widths = layout_widths
))

output_pdf <- file.path(fig_dir, "Figure08.pdf")
output_png <- file.path(fig_dir, "Figure08.png")

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

cropped_pdf <- tempfile(pattern = "Figure08_crop_", tmpdir = fig_dir, fileext = ".pdf")
crop_status <- suppressWarnings(system2("pdfcrop", args = c(output_pdf, cropped_pdf)))
if (identical(crop_status, 0L) && file.exists(cropped_pdf)) {
  ok_rename <- file.rename(cropped_pdf, output_pdf)
  if (!isTRUE(ok_rename)) {
    stop("pdfcrop succeeded but failed to replace the original PDF.")
  }
} else if (file.exists(cropped_pdf)) {
  unlink(cropped_pdf)
}

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
