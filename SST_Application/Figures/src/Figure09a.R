#!/usr/bin/env Rscript

# Figure 09a: posterior mean contour panels for the selected SST ensemble run.
# Left: posterior mean of X.
# Right: posterior mean of beta.

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
  library(grid)
  library(gridExtra)
})

source(file.path("modules", "plot_utils.R"))
source(file.path("modules", "plot_sst_maps.R"))
source(file.path(project_root, "Figures", "src", "figure_io_utils.R"))

figure_root <- file.path(project_root, "Figures")
best_rdata <- resolve_figure_chain_file(figure_root, "SST_FIGURE09A_CHAIN_FILE")
fig_dir <- resolve_figure_output_dir(figure_root, chain_file = best_rdata)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
grid_file <- file.path("Data", "processed", "pacific_sst_grid.csv")
land_file <- file.path("Data", "processed", "pacific_land_mask.csv")

for (f in c(best_rdata, grid_file, land_file)) {
  if (!file.exists(f)) {
    stop("Missing required file: ", f)
  }
}

result_env <- new.env(parent = emptyenv())
load(best_rdata, envir = result_env)
if (!exists("mcmc_artifact", envir = result_env, inherits = FALSE)) {
  stop("Best-run file does not contain mcmc_artifact: ", best_rdata)
}

artifact <- result_env$mcmc_artifact
output <- artifact$output
setup <- artifact$setup

required_output_fields <- c("X_post_mean", "beta_post_mean")
if (!all(required_output_fields %in% names(output))) {
  stop("Malformed mcmc_artifact$output: missing posterior mean fields.")
}

grid_df <- utils::read.csv(grid_file, stringsAsFactors = FALSE)
land_raw <- utils::read.csv(land_file, stringsAsFactors = FALSE)

required_grid_cols <- c("r", "c", "s", "lat", "lon")
if (!all(required_grid_cols %in% names(grid_df))) {
  stop("Grid file must contain columns: ", paste(required_grid_cols, collapse = ", "))
}
if (!all(c("s", "is_land") %in% names(land_raw))) {
  stop("Land-mask file must contain columns: s, is_land")
}

grid_df <- grid_df[order(grid_df$s), required_grid_cols]
land_raw <- land_raw[order(land_raw$s), c("s", "lat", "lon", "is_land")]

if (!identical(as.integer(grid_df$s), as.integer(land_raw$s))) {
  stop("Grid and land-mask files do not align by site id.")
}

R <- length(unique(grid_df$r))
C <- length(unique(grid_df$c))
S_full <- nrow(grid_df)
if (R * C != S_full) {
  stop("Grid dimensions inconsistent: R*C != number of rows in grid file.")
}

Xp <- matrix(grid_df$lon, nrow = R, ncol = C)
Yp <- matrix(grid_df$lat, nrow = R, ncol = C)
land_df <- data.frame(
  X = land_raw$lon,
  Y = land_raw$lat,
  is_land = as.logical(land_raw$is_land)
)
world_df <- build_world_boundaries(
  lon_limits = range(grid_df$lon, na.rm = TRUE),
  lat_limits = range(grid_df$lat, na.rm = TRUE)
)

ocean_sites_global <- as.integer(setup$ocean_sites_global)
if (length(ocean_sites_global) == 0L) {
  stop("Artifact setup is missing ocean_sites_global.")
}

expand_ocean_to_full <- function(vec_ocean, ocean_sites, S) {
  out <- rep(NA_real_, S)
  out[ocean_sites] <- as.numeric(vec_ocean)
  out
}

X_post_mean_full <- expand_ocean_to_full(output$X_post_mean, ocean_sites_global, S_full)
beta_post_mean_full <- expand_ocean_to_full(output$beta_post_mean, ocean_sites_global, S_full)

is_land <- as.logical(land_raw$is_land)
X_post_mean_full[is_land] <- NA_real_
beta_post_mean_full[is_land] <- NA_real_

X_post_mean_mat <- matrix(X_post_mean_full, nrow = R, ncol = C)
beta_post_mean_mat <- matrix(beta_post_mean_full, nrow = R, ncol = C)

palette_nondiv <- c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026")
palette_div_last_row <- c("#2166ac", "white", "#b2182b")
x_limits <- range(X_post_mean_full, na.rm = TRUE)
beta_limits <- c(-4, 4)
x_bins <- 45
beta_bins <- 100

panel_theme <- ggplot2::theme_bw(base_size = 9, base_family = "serif") +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = "black", linewidth = 0.35),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 8, color = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.25),
    plot.title = ggplot2::element_text(hjust = 0),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 8.5),
    legend.text = ggplot2::element_text(size = 7.5),
    legend.key.height = grid::unit(34, "pt"),
    legend.key.width = grid::unit(8, "pt"),
    legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    legend.box.spacing = grid::unit(6, "pt"),
    plot.margin = ggplot2::margin(3, 2, 3, 4, unit = "pt")
  )

x_plot <- suppressWarnings(plot_sst_contour_map(
  Z = X_post_mean_mat,
  Xp = Xp,
  Yp = Yp,
  bins = x_bins,
  palette = palette_nondiv,
  limits = x_limits,
  title = "(A) Posterior mean of X",
  land_df = land_df,
  world_df = world_df
)) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(0.12, "in"),
      barheight = grid::unit(2.0, "in")
    )
  ) +
  ggplot2::labs(fill = "deg C") +
  panel_theme

beta_plot <- suppressWarnings(plot_sst_contour_map(
  Z = beta_post_mean_mat,
  Xp = Xp,
  Yp = Yp,
  bins = beta_bins,
  palette = palette_div_last_row,
  zero_white = TRUE,
  limits = beta_limits,
  title = expression("(B) Posterior mean of " * beta),
  land_df = land_df,
  world_df = world_df
)) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(0.12, "in"),
      barheight = grid::unit(2.0, "in")
    )
  ) +
  ggplot2::labs(fill = "deg C") +
  panel_theme

figure_grob <- suppressWarnings(gridExtra::arrangeGrob(
  x_plot,
  beta_plot,
  nrow = 1,
  ncol = 2
))

output_pdf <- file.path(fig_dir, "Figure09a.pdf")
output_png <- file.path(fig_dir, "Figure09a.png")

figure_width_in <- 9.0
figure_height_in <- 3.9
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

cropped_pdf <- tempfile(pattern = "Figure09a_crop_", tmpdir = fig_dir, fileext = ".pdf")
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
