#!/usr/bin/env Rscript

resolve_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) return(NULL)
  normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE)
}

resolve_project_root <- function() {
  script_path <- resolve_script_path()
  candidates <- c(
    if (!is.null(script_path)) normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE) else NA_character_,
    normalizePath(getwd(), mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
  )
  candidates <- unique(candidates[!is.na(candidates)])

  for (cand in candidates) {
    if (file.exists(file.path(cand, "modules", "spectral2_ens.R")) &&
        dir.exists(file.path(cand, "Data"))) {
      return(cand)
    }
  }

  stop("Could not locate project root containing modules/spectral2_ens.R and Data/.", call. = FALSE)
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

resolve_api_file <- function(project_root, artifact) {
  data_dir_from_artifact <- artifact$setup$data_dir %||% ""
  input_dir_from_artifact <- artifact$setup$input_dir %||% ""
  candidate_base_dirs <- unique(Filter(nzchar, c(
    Sys.getenv("OXY_DATA_DIR", unset = ""),
    data_dir_from_artifact,
    file.path(project_root, "Data"),
    "Data"
  )))
  candidate_dirs <- unique(c(
    Sys.getenv("OXY_INPUT_DIR", unset = ""),
    input_dir_from_artifact,
    file.path(candidate_base_dirs, "inputs"),
    candidate_base_dirs
  ))
  candidate_dirs <- candidate_dirs[nzchar(candidate_dirs)]
  candidate_files <- unique(Filter(nzchar, c(
    Sys.getenv("OXY_API_FILE", unset = ""),
    file.path(candidate_dirs, "api.txt")
  )))
  hits <- candidate_files[file.exists(candidate_files)]
  if (length(hits) == 0L) {
    stop(
      "Missing api index file. Tried: ",
      paste(candidate_files, collapse = ", "),
      call. = FALSE
    )
  }

  normalizePath(hits[1L], winslash = "/", mustWork = FALSE)
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

`%||%` <- function(a, b) if (!is.null(a)) a else b

source(file.path(project_root, "Figures", "src", "plot_utils.R"))

extract_legend <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1L)) == "guide-box")
  if (length(idx) == 0L) return(NULL)
  g$grobs[[idx[1L]]]
}

figure_root <- file.path(project_root, "Figures")
chain_file <- resolve_figure_chain_file(
  figure_root,
  "OXY_FIGURE04B_CHAIN_FILE"
)
figure_dir <- resolve_figure_output_dir(figure_root, chain_file)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

ensemble_file <- Sys.getenv(
  "FK_ENSEMBLE_FILE",
  unset = file.path(project_root, "Data", "inputs", "oxy_ensemble_20_members.rds")
)

if (!file.exists(chain_file)) {
  stop("Missing chain file: ", chain_file, call. = FALSE)
}
if (!file.exists(ensemble_file)) {
  stop("Missing ensemble file: ", ensemble_file, call. = FALSE)
}

chain_env <- new.env(parent = emptyenv())
load(chain_file, envir = chain_env)
if (!exists("mcmc_artifact", envir = chain_env, inherits = FALSE)) {
  stop("Expected `mcmc_artifact` in chain file: ", chain_file, call. = FALSE)
}
artifact <- get("mcmc_artifact", envir = chain_env, inherits = FALSE)
api_file <- resolve_api_file(project_root, artifact)

required_artifact_fields <- c("output", "setup", "params")
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

X_post_sd <- compute_row_sds(X_chain[, keep_idx, drop = FALSE])
beta_post_sd <- compute_row_sds(beta_chain[, keep_idx, drop = FALSE])

ens <- readRDS(ensemble_file)
required_ensemble_fields <- c("domain", "members", "fields")
if (!all(required_ensemble_fields %in% names(ens))) {
  stop(
    "Ensemble RDS is missing required fields. Expected: ",
    paste(required_ensemble_fields, collapse = ", "),
    call. = FALSE
  )
}

lat_keep <- as.numeric(ens$domain$lat)
lon_keep <- as.numeric(ens$domain$lon)
R <- length(lat_keep)
C <- length(lon_keep)
S <- R * C
if (length(X_post_sd) != S || length(beta_post_sd) != S) {
  stop(
    sprintf(
      "Posterior SD lengths do not match grid size %d x %d = %d.",
      R, C, S
    ),
    call. = FALSE
  )
}

X_post_sd_mat <- matrix(X_post_sd, nrow = R, ncol = C)
beta_post_sd_mat <- matrix(beta_post_sd, nrow = R, ncol = C)
sd_diff_mat <- X_post_sd_mat - beta_post_sd_mat

obs_sites_global <- scan(api_file, what = integer(), quiet = TRUE)
if (length(obs_sites_global) == 0L) {
  stop("No observation indices found in api file: ", api_file, call. = FALSE)
}
if (any(!is.finite(obs_sites_global)) || any(obs_sites_global < 1L) || any(obs_sites_global > S)) {
  stop("Observation indices from api file are out of bounds for the OXY grid.", call. = FALSE)
}

grid_df <- data.frame(
  lon = rep(lon_keep, each = R),
  lat = rep(lat_keep, times = C)
)
obs_sites_df <- unique(grid_df[obs_sites_global, c("lon", "lat"), drop = FALSE])

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

shared_sd_limit <- max(c(X_post_sd_mat, beta_post_sd_mat), na.rm = TRUE)
if (!is.finite(shared_sd_limit) || shared_sd_limit <= 0) {
  shared_sd_limit <- 1
}
shared_sd_limit <- ceiling(shared_sd_limit * 10) / 10
shared_fill_limits <- c(0, shared_sd_limit)
shared_fill_breaks <- pretty(shared_fill_limits, n = 5)
shared_fill_breaks <- shared_fill_breaks[
  shared_fill_breaks >= shared_fill_limits[1] &
    shared_fill_breaks <= shared_fill_limits[2]
]

build_sd_fill_scale <- function() {
  ggplot2::scale_fill_gradientn(
    colors = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    limits = shared_fill_limits,
    breaks = shared_fill_breaks,
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
    barwidth = if (is_horizontal) {
      grid::unit(4.6, "in")
    } else {
      grid::unit(0.18, "in")
    },
    barheight = if (is_horizontal) {
      grid::unit(0.18, "in")
    } else {
      grid::unit(1.85, "in")
    },
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
      bins = 45
    ) +
    ggplot2::coord_equal(xlim = plot_lon_limits, ylim = plot_lat_limits, expand = FALSE) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = plot_lon_limits) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = plot_lat_limits) +
    build_sd_fill_scale() +
    ggplot2::labs(title = title, fill = NULL) +
    build_oxy_panel_theme() +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) legend_position else "none",
      legend.direction = legend_direction,
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    )

  if (isTRUE(show_legend)) {
    panel <- panel + ggplot2::guides(fill = build_sd_fill_guide(legend_direction))
  }

  panel <- add_world_boundaries(panel, world_df = world_df)
  if (nrow(obs_sites_df) > 0L) {
    panel <- panel +
      ggplot2::geom_point(
        data = obs_sites_df,
        ggplot2::aes(x = lon, y = lat),
        inherit.aes = FALSE,
        shape = 16,
        size = 0.5,
        color = "black",
        alpha = 0.9
      )
  }

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
    build_oxy_panel_theme() +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) "right" else "none",
      legend.direction = "vertical"
    )

  if (isTRUE(show_legend)) {
    panel <- panel + ggplot2::guides(fill = build_diff_fill_guide())
  }

  panel <- add_world_boundaries(panel, world_df = world_df)
  if (nrow(obs_sites_df) > 0L) {
    panel <- panel +
      ggplot2::geom_point(
        data = obs_sites_df,
        ggplot2::aes(x = lon, y = lat),
        inherit.aes = FALSE,
        shape = 16,
        size = 0.5,
        color = "black",
        alpha = 0.9
      )
  }

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

probe_pdf <- tempfile(fileext = ".pdf")
grDevices::pdf(file = probe_pdf, width = 13.5, height = 3.8)
sd_legend_grob <- suppressWarnings(extract_legend(build_sd_panel(
  z_mat = X_post_sd_mat,
  title = NULL,
  show_legend = TRUE,
  legend_position = "left",
  legend_direction = "vertical"
)))
if (is.null(sd_legend_grob)) {
  stop("Could not extract common SD legend for Figure04b.", call. = FALSE)
}

diff_legend_grob <- suppressWarnings(extract_legend(build_diff_panel(
  z_mat = sd_diff_mat,
  title = NULL,
  show_legend = TRUE
)))
if (is.null(diff_legend_grob)) {
  stop("Could not extract SD-difference legend for Figure04b.", call. = FALSE)
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
grDevices::dev.off()
unlink(probe_pdf)

output_pdf <- file.path(figure_dir, "Figure04b.pdf")
output_png <- file.path(figure_dir, "Figure04b.png")

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
