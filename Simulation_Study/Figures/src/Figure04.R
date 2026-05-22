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

figure_root <- file.path(project_root, "Figures")
chain_file <- resolve_figure_chain_file(
  figure_root,
  "OXY_FIGURE04_CHAIN_FILE"
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

required_artifact_fields <- c("output", "setup", "params")
if (!all(required_artifact_fields %in% names(artifact))) {
  stop(
    "Chain artifact is missing required fields. Expected: ",
    paste(required_artifact_fields, collapse = ", "),
    call. = FALSE
  )
}

X_post_mean <- as.numeric(artifact$output$X_post_mean %||% numeric(0))
beta_post_mean <- as.numeric(artifact$output$beta_post_mean %||% numeric(0))
if (length(X_post_mean) == 0L || length(beta_post_mean) == 0L) {
  stop("Posterior means are missing from chain artifact: ", chain_file, call. = FALSE)
}

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
if (length(X_post_mean) != S || length(beta_post_mean) != S) {
  stop(
    sprintf(
      "Posterior mean lengths do not match grid size %d x %d = %d.",
      R, C, S
    ),
    call. = FALSE
  )
}

X_post_mat <- matrix(X_post_mean, nrow = R, ncol = C)
beta_post_mat <- matrix(beta_post_mean, nrow = R, ncol = C)

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

x_fill_limits <- range(c(ens$fields, X_post_mat), na.rm = TRUE)
beta_limit <- max(abs(beta_post_mat), na.rm = TRUE)
beta_limit <- ceiling(beta_limit * 10) / 10

x_plot <- plot_oxy_contour(
  Z = X_post_mat,
  Xp = Xp,
  Yp = Yp,
  title = expression("(A) Posterior mean of " * X),
  fill_label = NULL,
  land_df = NULL,
  world_df = world_df,
  x_limits = plot_lon_limits,
  y_limits = plot_lat_limits,
  fill_limits = x_fill_limits
) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      direction = "vertical",
      title = NULL,
      label.position = "right",
      barwidth = grid::unit(0.18, "in"),
      barheight = grid::unit(2.6, "in"),
      ticks.colour = "black",
      frame.colour = "black",
      label.theme = ggplot2::element_text(size = 8, color = "black")
    )
  ) +
  ggplot2::theme(
    legend.position = "right",
    legend.direction = "vertical"
  )

beta_df <- data.frame(
  X = as.vector(Xp),
  Y = as.vector(Yp),
  Z = as.vector(beta_post_mat)
)

beta_plot <- ggplot2::ggplot(beta_df, ggplot2::aes(x = X, y = Y, z = Z)) +
  ggplot2::geom_contour_filled(
    ggplot2::aes(fill = after_stat((level_low + level_high) / 2)),
    bins = 45
  ) +
  ggplot2::coord_equal(xlim = plot_lon_limits, ylim = plot_lat_limits, expand = FALSE) +
  ggplot2::scale_x_continuous(expand = c(0, 0), limits = plot_lon_limits) +
  ggplot2::scale_y_continuous(expand = c(0, 0), limits = plot_lat_limits) +
  ggplot2::scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-beta_limit, beta_limit),
    oob = scales::squish
  ) +
  ggplot2::labs(
    title = expression("(B) Posterior mean of " * beta),
    fill = NULL
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      direction = "vertical",
      title = NULL,
      label.position = "right",
      barwidth = grid::unit(0.18, "in"),
      barheight = grid::unit(2.6, "in"),
      ticks.colour = "black",
      frame.colour = "black",
      label.theme = ggplot2::element_text(size = 8, color = "black")
    )
  ) +
  build_oxy_panel_theme() +
  ggplot2::theme(
    legend.position = "right",
    legend.direction = "vertical"
  )

beta_plot <- add_world_boundaries(beta_plot, world_df = world_df)

probe_pdf <- tempfile(fileext = ".pdf")
grDevices::pdf(file = probe_pdf, width = 7, height = 7)
panel_grobs <- lapply(list(x_plot, beta_plot), ggplot2::ggplotGrob)
common_widths <- do.call(grid::unit.pmax, lapply(panel_grobs, function(g) g$widths))
common_heights <- do.call(grid::unit.pmax, lapply(panel_grobs, function(g) g$heights))
panel_grobs <- lapply(panel_grobs, function(g) {
  g$widths <- common_widths
  g$heights <- common_heights
  g
})

figure_grob <- suppressWarnings(gridExtra::arrangeGrob(
  grobs = panel_grobs,
  nrow = 1,
  ncol = 2,
  padding = grid::unit(0, "pt")
))
grDevices::dev.off()
unlink(probe_pdf)

output_pdf <- file.path(figure_dir, "Figure04.pdf")
output_png <- file.path(figure_dir, "Figure04.png")

pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
suppressWarnings(ggplot2::ggsave(
  filename = output_pdf,
  plot = figure_grob,
  device = pdf_device,
  width = 10.2,
  height = 4.8,
  units = "in",
  bg = "white"
))
crop_pdf_in_place(output_pdf)

suppressWarnings(ggplot2::ggsave(
  filename = output_png,
  plot = figure_grob,
  width = 10.2,
  height = 4.8,
  units = "in",
  dpi = 600,
  bg = "white"
))

message("Chain file: ", normalizePath(chain_file, mustWork = FALSE))
message("Saved figure PDF: ", normalizePath(output_pdf, mustWork = FALSE))
message("Saved figure PNG: ", normalizePath(output_png, mustWork = FALSE))
