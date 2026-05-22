#!/usr/bin/env Rscript

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
  library(grid)
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
true_file <- Sys.getenv(
  "OXY_TRUE_FILE",
  unset = file.path(input_dir, "trueXX.txt")
)
api_file <- Sys.getenv(
  "OXY_API_FILE",
  unset = file.path(input_dir, "api.txt")
)

if (!file.exists(ensemble_file)) {
  stop("Missing ensemble file: ", ensemble_file, call. = FALSE)
}
if (!file.exists(true_file)) {
  stop("Missing true field file: ", true_file, call. = FALSE)
}
if (!file.exists(api_file)) {
  stop("Missing api index file: ", api_file, call. = FALSE)
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

true_mat <- as.matrix(utils::read.table(true_file))
if (!identical(dim(true_mat), c(R, C))) {
  stop(
    sprintf(
      "Expected true field dimensions %d x %d, found %d x %d.",
      R, C, nrow(true_mat), ncol(true_mat)
    ),
    call. = FALSE
  )
}

api_idx <- scan(api_file, what = integer(), quiet = TRUE)
if (length(api_idx) == 0L) {
  stop("No observation indices found in api file.", call. = FALSE)
}
if (any(!is.finite(api_idx)) || any(api_idx < 1L) || any(api_idx > R * C)) {
  stop("Observation indices must be integers between 1 and 703.", call. = FALSE)
}

mean_mat <- apply(ens$fields, c(1, 2), mean, na.rm = TRUE)
bias_mat <- mean_mat - true_mat

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

obs_df <- data.frame(
  X = as.vector(Xp)[api_idx],
  Y = as.vector(Yp)[api_idx]
)

true_title <- expression("(A) " * X)
bias_title <- expression("(B) " * bar(Y)^c - X)

true_plot <- plot_oxy_contour(
  Z = true_mat,
  Xp = Xp,
  Yp = Yp,
  title = true_title,
  fill_label = NULL,
  land_df = NULL,
  world_df = world_df,
  x_limits = plot_lon_limits,
  y_limits = plot_lat_limits,
  fill_limits = range(true_mat, na.rm = TRUE)
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
  ) +
  ggplot2::geom_point(
    data = obs_df,
    mapping = ggplot2::aes(x = X, y = Y),
    inherit.aes = FALSE,
    color = "black",
    shape = 16,
    size = 0.42,
    alpha = 0.8
  )

bias_df <- data.frame(
  X = as.vector(Xp),
  Y = as.vector(Yp),
  Z = as.vector(bias_mat)
)

bias_plot <- ggplot2::ggplot(bias_df, ggplot2::aes(x = X, y = Y, z = Z)) +
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
    limits = c(-10, 10),
    oob = scales::squish
  ) +
  ggplot2::labs(
    title = bias_title,
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

bias_plot <- add_world_boundaries(bias_plot, world_df = world_df)

panel_grobs <- lapply(list(true_plot, bias_plot), ggplot2::ggplotGrob)
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

output_dir <- Sys.getenv("FK_OUTPUT_DIR", unset = file.path(figure_dir, "out"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_pdf <- file.path(output_dir, "Figure03.pdf")
output_png <- file.path(output_dir, "Figure03.png")

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

message("Saved figure PDF: ", normalizePath(output_pdf, mustWork = FALSE))
message("Saved figure PNG: ", normalizePath(output_png, mustWork = FALSE))
