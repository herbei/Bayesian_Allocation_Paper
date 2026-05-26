#!/usr/bin/env Rscript

# Figure 09c: posterior SD panels for X and beta, plus their sitewise difference.

resolve_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) return(NULL)
  normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE)
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

  stop("Could not locate project root containing modules/plot_sst_maps.R", call. = FALSE)
}

compute_row_sds <- function(draw_matrix) {
  draw_matrix <- as.matrix(draw_matrix)
  n_draws <- ncol(draw_matrix)
  if (n_draws == 0L) {
    stop("Cannot compute posterior standard deviations with zero draws.", call. = FALSE)
  }
  if (n_draws == 1L) {
    return(rep(0, nrow(draw_matrix)))
  }

  draw_means <- rowMeans(draw_matrix)
  centered_draws <- draw_matrix - draw_means
  sqrt(pmax(0, rowSums(centered_draws^2) / (n_draws - 1L)))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

resolve_keep_idx <- function(artifact, store) {
  keep_idx <- as.integer(artifact$output$keep_idx %||% integer(0))
  n_store <- ncol(as.matrix(store$X_chain %||% matrix(numeric(0), nrow = 0L, ncol = 0L)))

  keep_idx <- keep_idx[
    is.finite(keep_idx) &
      keep_idx >= 1L &
      keep_idx <= n_store
  ]
  if (length(keep_idx) > 0L) {
    return(unique(keep_idx))
  }

  iter_chain <- as.integer(store$iter_chain %||% integer(0))
  burn_in <- as.integer(artifact$run_config$burn_in %||% 0L)
  if (length(iter_chain) == n_store) {
    keep_idx <- which(iter_chain > burn_in)
    if (length(keep_idx) > 0L) {
      return(keep_idx)
    }
  }

  seq_len(n_store)
}

expand_ocean_to_full <- function(vec_ocean, ocean_sites, S_full) {
  out <- rep(NA_real_, S_full)
  out[ocean_sites] <- as.numeric(vec_ocean)
  out
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

extract_legend <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1L)) == "guide-box")
  if (length(idx) == 0L) return(NULL)
  g$grobs[[idx[1L]]]
}

build_sst_panel_theme <- function() {
  ggplot2::theme_bw(base_size = 9, base_family = "serif") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.35),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 8, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25),
      plot.title = ggplot2::element_text(hjust = 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
      plot.margin = ggplot2::margin(3, 3, 3, 3, unit = "pt")
    )
}

figure_root <- file.path(project_root, "Figures")
chain_file <- resolve_figure_chain_file(
  figure_root,
  c("SST_FIGURE09C_SD_CHAIN_FILE", "SST_FIGURE09C_CHAIN_FILE")
)
figure_dir <- resolve_figure_output_dir(figure_root, chain_file)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

grid_file <- file.path("Data", "processed", "pacific_sst_grid.csv")
land_file <- file.path("Data", "processed", "pacific_land_mask.csv")

for (f in c(chain_file, grid_file, land_file)) {
  if (!file.exists(f)) {
    stop("Missing required file: ", f, call. = FALSE)
  }
}

chain_env <- new.env(parent = emptyenv())
load(chain_file, envir = chain_env)
if (!exists("mcmc_artifact", envir = chain_env, inherits = FALSE)) {
  stop("Expected `mcmc_artifact` in chain file: ", chain_file, call. = FALSE)
}
artifact <- get("mcmc_artifact", envir = chain_env, inherits = FALSE)

required_artifact_fields <- c("output", "setup")
if (!all(required_artifact_fields %in% names(artifact))) {
  stop(
    "Chain artifact is missing required fields. Expected: ",
    paste(required_artifact_fields, collapse = ", "),
    call. = FALSE
  )
}

store <- artifact$output$store %||% list()
required_store_fields <- c("X_chain", "beta_chain")
if (!all(required_store_fields %in% names(store))) {
  stop(
    "Chain artifact store is missing required fields. Expected: ",
    paste(required_store_fields, collapse = ", "),
    call. = FALSE
  )
}

X_chain <- as.matrix(store$X_chain)
beta_chain <- as.matrix(store$beta_chain)
if (ncol(X_chain) == 0L || ncol(beta_chain) == 0L) {
  stop("Chain artifact does not contain stored posterior draws.", call. = FALSE)
}
if (!identical(dim(X_chain), dim(beta_chain))) {
  stop("Stored X and beta draw matrices do not have matching dimensions.", call. = FALSE)
}

keep_idx <- resolve_keep_idx(artifact, store)
if (length(keep_idx) == 0L) {
  stop("No kept posterior draws are available in the chain artifact.", call. = FALSE)
}

grid_df <- utils::read.csv(grid_file, stringsAsFactors = FALSE)
land_raw <- utils::read.csv(land_file, stringsAsFactors = FALSE)

required_grid_cols <- c("r", "c", "s", "lat", "lon")
if (!all(required_grid_cols %in% names(grid_df))) {
  stop("Grid file must contain columns: ", paste(required_grid_cols, collapse = ", "), call. = FALSE)
}
if (!all(c("s", "lat", "lon", "is_land") %in% names(land_raw))) {
  stop("Land-mask file must contain columns: s, lat, lon, is_land", call. = FALSE)
}

grid_df <- grid_df[order(grid_df$s), required_grid_cols]
land_raw <- land_raw[order(land_raw$s), c("s", "lat", "lon", "is_land")]

if (!identical(as.integer(grid_df$s), as.integer(land_raw$s))) {
  stop("Grid and land-mask files do not align by site id.", call. = FALSE)
}

R <- length(unique(grid_df$r))
C <- length(unique(grid_df$c))
S_full <- nrow(grid_df)
if (R * C != S_full) {
  stop("Grid dimensions inconsistent: R*C != number of rows in grid file.", call. = FALSE)
}

ocean_sites_global <- as.integer(artifact$setup$ocean_sites_global %||% integer(0))
if (length(ocean_sites_global) == 0L) {
  stop("Artifact setup is missing ocean_sites_global.", call. = FALSE)
}
if (nrow(X_chain) != length(ocean_sites_global) || nrow(beta_chain) != length(ocean_sites_global)) {
  stop("Posterior chain rows do not align with setup$ocean_sites_global.", call. = FALSE)
}

Xp <- matrix(grid_df$lon, nrow = R, ncol = C)
Yp <- matrix(grid_df$lat, nrow = R, ncol = C)
land_df <- data.frame(
  X = land_raw$lon,
  Y = land_raw$lat,
  is_land = as.logical(land_raw$is_land)
)
plot_lon_limits <- range(grid_df$lon, na.rm = TRUE)
plot_lat_limits <- range(grid_df$lat, na.rm = TRUE)
world_df <- build_world_boundaries(
  lon_limits = plot_lon_limits,
  lat_limits = plot_lat_limits,
  lon_buffer = 0,
  lat_buffer = 0
)

X_post_sd <- compute_row_sds(X_chain[, keep_idx, drop = FALSE])
beta_post_sd <- compute_row_sds(beta_chain[, keep_idx, drop = FALSE])

X_post_sd_full <- expand_ocean_to_full(X_post_sd, ocean_sites_global, S_full)
beta_post_sd_full <- expand_ocean_to_full(beta_post_sd, ocean_sites_global, S_full)

is_land <- as.logical(land_raw$is_land)
X_post_sd_full[is_land] <- NA_real_
beta_post_sd_full[is_land] <- NA_real_

X_post_sd_mat <- matrix(X_post_sd_full, nrow = R, ncol = C)
beta_post_sd_mat <- matrix(beta_post_sd_full, nrow = R, ncol = C)
sd_diff_mat <- X_post_sd_mat - beta_post_sd_mat

shared_fill_limits <- c(0, 1)
shared_fill_breaks <- seq(0, 1, by = 0.2)
shared_fill_labels <- c("0", "0.2", "0.4", "0.6", "0.8", ">1")
sd_bins <- 220L

build_sd_fill_scale <- function() {
  ggplot2::scale_fill_gradientn(
    colors = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    limits = shared_fill_limits,
    breaks = shared_fill_breaks,
    labels = shared_fill_labels,
    oob = scales::squish
  )
}

build_sd_fill_guide <- function(direction = "vertical") {
  is_horizontal <- identical(direction, "horizontal")

  ggplot2::guide_colorbar(
    direction = direction,
    title = NULL,
    title.position = "top",
    label.position = if (is_horizontal) "bottom" else "left",
    barwidth = if (is_horizontal) grid::unit(4.6, "in") else grid::unit(0.18, "in"),
    barheight = if (is_horizontal) grid::unit(0.18, "in") else grid::unit(1.85, "in"),
    ticks.colour = "black",
    frame.colour = "black",
    title.theme = ggplot2::element_text(size = 8, color = "black"),
    label.theme = ggplot2::element_text(size = 8, color = "black")
  )
}

build_sd_panel <- function(
  z_mat,
  title,
  show_legend = TRUE,
  legend_position = "right",
  legend_direction = "vertical"
) {
  plot_df <- data.frame(
    X = as.vector(Xp),
    Y = as.vector(Yp),
    Z = as.vector(z_mat)
  )

  panel <- ggplot2::ggplot(plot_df, ggplot2::aes(x = X, y = Y, z = Z)) +
    ggplot2::geom_contour_filled(
      ggplot2::aes(fill = after_stat((level_low + level_high) / 2)),
      bins = sd_bins
    ) +
    ggplot2::coord_equal(xlim = plot_lon_limits, ylim = plot_lat_limits, expand = FALSE) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = plot_lon_limits) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = plot_lat_limits) +
    build_sd_fill_scale() +
    ggplot2::labs(title = title, fill = NULL) +
    build_sst_panel_theme() +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) legend_position else "none",
      legend.direction = legend_direction
    )

  if (isTRUE(show_legend)) {
    panel <- panel + ggplot2::guides(fill = build_sd_fill_guide(legend_direction))
  }

  panel <- add_land_mask(panel, land_df = land_df, fill = "grey80")
  panel <- add_world_boundaries(panel, world_df = world_df)

  panel
}

diff_limit <- max(abs(sd_diff_mat), na.rm = TRUE)
if (!is.finite(diff_limit) || diff_limit <= 0) {
  diff_limit <- 0.01
}
diff_limit <- max(0.01, ceiling(diff_limit * 100) / 100)
diff_fill_limits <- c(-diff_limit, diff_limit)
diff_fill_breaks <- seq(diff_fill_limits[1], diff_fill_limits[2], length.out = 7L)

build_diff_fill_scale <- function() {
  ggplot2::scale_fill_gradientn(
    colors = c("#8c1d40", "#d9899b", "#f7f7f7", "#8db9da", "#1d4e89"),
    limits = diff_fill_limits,
    breaks = diff_fill_breaks,
    labels = formatC(diff_fill_breaks, format = "f", digits = 2),
    oob = scales::squish
  )
}

build_diff_fill_guide <- function() {
  ggplot2::guide_colorbar(
    direction = "vertical",
    title = NULL,
    title.position = "top",
    label.position = "right",
    barwidth = grid::unit(0.18, "in"),
    barheight = grid::unit(1.85, "in"),
    ticks.colour = "black",
    frame.colour = "black",
    title.theme = ggplot2::element_text(size = 8, color = "black"),
    label.theme = ggplot2::element_text(size = 8, color = "black")
  )
}

build_diff_panel <- function(z_mat, title, show_legend = TRUE) {
  plot_df <- data.frame(
    X = as.vector(Xp),
    Y = as.vector(Yp),
    Z = as.vector(z_mat)
  )

  panel <- ggplot2::ggplot(plot_df, ggplot2::aes(x = X, y = Y, z = Z)) +
    ggplot2::geom_contour_filled(
      ggplot2::aes(fill = after_stat((level_low + level_high) / 2)),
      bins = 40
    ) +
    ggplot2::geom_contour(
      breaks = 0,
      color = "black",
      linewidth = 0.35,
      alpha = 0.8
    ) +
    ggplot2::coord_equal(xlim = plot_lon_limits, ylim = plot_lat_limits, expand = FALSE) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = plot_lon_limits) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = plot_lat_limits) +
    build_diff_fill_scale() +
    ggplot2::labs(title = title, fill = NULL) +
    build_sst_panel_theme() +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) "right" else "none",
      legend.direction = "vertical"
    )

  if (isTRUE(show_legend)) {
    panel <- panel + ggplot2::guides(fill = build_diff_fill_guide())
  }

  panel <- add_land_mask(panel, land_df = land_df, fill = "grey80")
  panel <- add_world_boundaries(panel, world_df = world_df)

  panel
}

x_sd_plot <- build_sd_panel(
  z_mat = X_post_sd_mat,
  title = expression("(A) Posterior SD of " * X),
  show_legend = FALSE
)

beta_sd_plot <- build_sd_panel(
  z_mat = beta_post_sd_mat,
  title = expression("(B) Posterior SD of " * beta),
  show_legend = FALSE
)

diff_sd_plot <- build_diff_panel(
  z_mat = sd_diff_mat,
  title = expression("(C) Posterior SD difference: SD(" * X * ") - SD(" * beta * ")"),
  show_legend = FALSE
)

sd_legend_grob <- suppressWarnings(extract_legend(build_sd_panel(
  z_mat = X_post_sd_mat,
  title = NULL,
  show_legend = TRUE,
  legend_position = "left",
  legend_direction = "vertical"
)))
if (is.null(sd_legend_grob)) {
  stop("Could not extract common SD legend.", call. = FALSE)
}

diff_legend_grob <- suppressWarnings(extract_legend(build_diff_panel(
  z_mat = sd_diff_mat,
  title = NULL,
  show_legend = TRUE
)))
if (is.null(diff_legend_grob)) {
  stop("Could not extract SD-difference legend.", call. = FALSE)
}

panel_grobs <- lapply(list(x_sd_plot, beta_sd_plot, diff_sd_plot), ggplot2::ggplotGrob)
common_widths <- do.call(grid::unit.pmax, lapply(panel_grobs, function(g) g$widths))
common_heights <- do.call(grid::unit.pmax, lapply(panel_grobs, function(g) g$heights))
panel_grobs <- lapply(panel_grobs, function(g) {
  g$widths <- common_widths
  g$heights <- common_heights
  g
})

c_panel_grob <- suppressWarnings(gridExtra::arrangeGrob(
  panel_grobs[[3]],
  grid::nullGrob(),
  diff_legend_grob,
  ncol = 3,
  widths = c(1, 0.04, 0.18),
  padding = grid::unit(0, "pt")
))

ab_panel_grob <- suppressWarnings(gridExtra::arrangeGrob(
  grobs = panel_grobs[1:2],
  nrow = 1,
  ncol = 2,
  padding = grid::unit(0, "pt")
))

left_block_grob <- suppressWarnings(gridExtra::arrangeGrob(
  sd_legend_grob,
  ab_panel_grob,
  ncol = 2,
  widths = c(0.09, 1),
  padding = grid::unit(0, "pt")
))

figure_grob <- suppressWarnings(gridExtra::arrangeGrob(
  left_block_grob,
  c_panel_grob,
  ncol = 2,
  widths = c(2.12, 1.18),
  padding = grid::unit(0, "pt")
))

output_pdf <- file.path(figure_dir, "Figure09c.pdf")
output_png <- file.path(figure_dir, "Figure09c.png")

pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
suppressWarnings(ggplot2::ggsave(
  filename = output_pdf,
  plot = figure_grob,
  device = pdf_device,
  width = 13.6,
  height = 4.15,
  units = "in",
  bg = "white"
))
crop_pdf_in_place(output_pdf)

suppressWarnings(ggplot2::ggsave(
  filename = output_png,
  plot = figure_grob,
  width = 13.6,
  height = 4.15,
  units = "in",
  dpi = 600,
  bg = "white"
))

message("Chain file: ", normalizePath(chain_file, mustWork = FALSE))
message("Saved figure PDF: ", normalizePath(output_pdf, mustWork = FALSE))
message("Saved figure PNG: ", normalizePath(output_png, mustWork = FALSE))
