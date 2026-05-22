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
figure_seed_value <- Sys.getenv("FK_FIGURE_SEED", unset = "1")
figure_seed <- if (nzchar(figure_seed_value)) as.integer(figure_seed_value) else NULL

if (!file.exists(ensemble_file)) {
  stop("Missing ensemble file: ", ensemble_file, call. = FALSE)
}

ensemble <- readRDS(ensemble_file)
required_fields <- c("meta", "domain", "members", "fields")
if (!all(required_fields %in% names(ensemble))) {
  stop(
    "Ensemble file is missing required fields. Expected: ",
    paste(required_fields, collapse = ", "),
    call. = FALSE
  )
}

n_members <- length(ensemble$members)
if (n_members < 4L) {
  stop("Need at least 4 ensemble members to build the panel figure.", call. = FALSE)
}

if (!is.null(figure_seed)) {
  set.seed(figure_seed)
}
selected_index <- sort(sample.int(n_members, size = 4L, replace = FALSE))
selected_members <- ensemble$members[selected_index]

grid_lon <- as.numeric(ensemble$domain$lon)
grid_lat <- as.numeric(ensemble$domain$lat)
n_rows <- length(grid_lat)
n_cols <- length(grid_lon)
Xp <- matrix(rep(grid_lon, each = n_rows), nrow = n_rows, ncol = n_cols)
Yp <- matrix(rep(grid_lat, times = n_cols), nrow = n_rows, ncol = n_cols)

plot_lon_limits <- if (!is.null(ensemble$domain$plot_lon_limits)) {
  as.numeric(ensemble$domain$plot_lon_limits)
} else {
  range(grid_lon, na.rm = TRUE) + c(-2, 2)
}
plot_lat_limits <- if (!is.null(ensemble$domain$plot_lat_limits)) {
  as.numeric(ensemble$domain$plot_lat_limits)
} else {
  range(grid_lat, na.rm = TRUE) + c(-2, 2)
}

world_df <- build_world_boundaries(
  lon_limits = plot_lon_limits,
  lat_limits = plot_lat_limits,
  lon_buffer = 0,
  lat_buffer = 0
)

common_limits <- range(ensemble$fields, na.rm = TRUE)
sst_panel_title_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(
    family = "serif",
    face = "plain",
    color = "black",
    size = 10.8,
    hjust = 0,
    margin = ggplot2::margin(0, 0, 4.5, 0, unit = "pt")
  )
)
panel_plots <- vector("list", length = 4L)
for (i in seq_along(selected_index)) {
  member_idx <- selected_index[i]
  panel <- plot_oxy_contour(
    Z = ensemble$fields[, , member_idx],
    Xp = Xp,
    Yp = Yp,
    title = sprintf("(%s) Member %s", LETTERS[i], ensemble$members[member_idx]),
    fill_label = NULL,
    land_df = NULL,
    world_df = world_df,
    x_limits = plot_lon_limits,
    y_limits = plot_lat_limits,
    fill_limits = common_limits
  ) + sst_panel_title_theme

  if (i <= 2L) {
    panel <- panel + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(1, 3, -5, 3, unit = "pt")
    )
  } else {
    panel <- panel + ggplot2::theme(
      plot.margin = ggplot2::margin(0, 3, 3, 3, unit = "pt")
    )
  }
  if (i %% 2L == 0L) {
    panel <- panel + ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  }
  panel_plots[[i]] <- panel
}

output_dir <- Sys.getenv("FK_OUTPUT_DIR", unset = file.path(figure_dir, "out"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_pdf <- file.path(output_dir, "Figure01.pdf")
output_png <- file.path(output_dir, "Figure01.png")
save_panel_figure(panel_plots = panel_plots, output_pdf = output_pdf, output_png = output_png)

message("Selected members: ", paste(selected_members, collapse = ", "))
message("Saved figure PDF: ", output_pdf)
message("Saved figure PNG: ", output_png)
