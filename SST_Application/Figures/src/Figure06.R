#!/usr/bin/env Rscript

# Figure 06: two-panel map figure.
# Left: processed observational SST contour with sampled-site locations.
# Right: full Pacific analysis grid on map.

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

grid_file <- file.path("Data", "processed", "pacific_sst_grid.csv")
obs_file <- file.path("Data", "processed", "pacific_obs_sites.csv")
land_file <- file.path("Data", "processed", "pacific_land_mask.csv")
edges_file <- file.path("Data", "processed", "pacific_ocean_edges.csv")
yo_full_file <- file.path("Data", "processed", "pacific_Yo_full.csv")
bundle_file <- file.path("Data", "processed", "data.rds")

for (f in c(grid_file, obs_file, land_file, edges_file, yo_full_file)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}

grid_df <- utils::read.csv(grid_file, stringsAsFactors = FALSE)
obs_df <- utils::read.csv(obs_file, stringsAsFactors = FALSE)
land_raw <- utils::read.csv(land_file, stringsAsFactors = FALSE)
edge_df <- utils::read.csv(edges_file, stringsAsFactors = FALSE)
yo_full_df <- utils::read.csv(yo_full_file, stringsAsFactors = FALSE)

required_grid_cols <- c("r", "c", "s", "lat", "lon")
if (!all(required_grid_cols %in% names(grid_df))) {
  stop("Grid file must contain columns: ", paste(required_grid_cols, collapse = ", "))
}
if (!all(c("s") %in% names(obs_df))) {
  stop("Obs-sites file must contain column: s")
}
if (!all(c("s", "lat", "lon", "is_land") %in% names(land_raw))) {
  stop("Land-mask file must contain columns: s, lat, lon, is_land")
}
if (!all(c("s_from", "s_to") %in% names(edge_df))) {
  stop("Ocean-edges file must contain columns: s_from, s_to")
}
if (!all(c("s", "Yo_full") %in% names(yo_full_df))) {
  stop("Full observation file must contain columns: s, Yo_full")
}

grid_df <- grid_df[order(grid_df$s), required_grid_cols]
land_raw <- land_raw[order(land_raw$s), c("s", "lat", "lon", "is_land")]
obs_df <- obs_df[order(obs_df$s), , drop = FALSE]
yo_full_df <- yo_full_df[order(yo_full_df$s), c("s", "Yo_full"), drop = FALSE]

if (!identical(as.integer(grid_df$s), as.integer(land_raw$s))) {
  stop("pacific_sst_grid.csv and pacific_land_mask.csv do not align by site id s.")
}
if (!identical(as.integer(grid_df$s), as.integer(yo_full_df$s))) {
  stop("pacific_sst_grid.csv and pacific_Yo_full.csv do not align by site id s.")
}

R <- length(unique(grid_df$r))
C <- length(unique(grid_df$c))
if (R * C != nrow(grid_df)) {
  stop("Grid dimensions inconsistent: R*C != number of rows in grid file.")
}

if (any(!obs_df$s %in% grid_df$s)) {
  stop("pacific_obs_sites.csv contains site ids not present in pacific_sst_grid.csv.")
}

Xp <- matrix(grid_df$lon, nrow = R, ncol = C)
Yp <- matrix(grid_df$lat, nrow = R, ncol = C)
land_df <- data.frame(X = land_raw$lon, Y = land_raw$lat, is_land = as.logical(land_raw$is_land))
world_df <- build_world_boundaries(
  lon_limits = range(grid_df$lon, na.rm = TRUE),
  lat_limits = range(grid_df$lat, na.rm = TRUE)
)

target_ym <- "2012-01"
if (file.exists(bundle_file)) {
  sst_bundle <- readRDS(bundle_file)
  if (!is.null(sst_bundle$meta$target_ym) &&
      nzchar(as.character(sst_bundle$meta$target_ym))) {
    target_ym <- as.character(sst_bundle$meta$target_ym)
  }
}

get_ym_oisst <- function(nc) {
  time <- ncdf4::ncvar_get(nc, "time")
  dates <- as.Date("1800-01-01") + time
  format(dates, "%Y-%m")
}

extract_window <- function(nc, varname, t_index, lat_keep, lon_keep, missing_value = NA_real_,
                           apply_scale = TRUE) {
  lat <- ncdf4::ncvar_get(nc, "lat")
  lon <- ncdf4::ncvar_get(nc, "lon")

  lat_idx <- match(lat_keep, lat)
  lon_idx <- match(lon_keep, lon)
  if (any(is.na(lat_idx))) stop("Some grid lat values were not found in OISST lat.")
  if (any(is.na(lon_idx))) stop("Some grid lon values were not found in OISST lon.")

  v <- nc$var[[varname]]
  if (is.null(v)) stop("Variable not found in OISST file: ", varname)

  dim_names <- vapply(v$dim, function(d) d$name, character(1))
  i_lon <- match("lon", dim_names)
  i_lat <- match("lat", dim_names)
  i_time <- match("time", dim_names)
  if (any(is.na(c(i_lon, i_lat, i_time)))) {
    stop("Expected lon/lat/time dimensions for variable ", varname)
  }

  start <- rep(1L, length(dim_names))
  count <- rep(1L, length(dim_names))
  start[i_lon] <- min(lon_idx)
  count[i_lon] <- length(lon_idx)
  start[i_lat] <- min(lat_idx)
  count[i_lat] <- length(lat_idx)
  start[i_time] <- t_index
  count[i_time] <- 1L

  x <- ncdf4::ncvar_get(
    nc,
    varname,
    start = start,
    count = count,
    collapse_degen = FALSE,
    raw_datavals = !apply_scale
  )

  x <- aperm(x, perm = c(i_lat, i_lon, setdiff(seq_along(dim_names), c(i_lat, i_lon))))
  x <- drop(x)
  if (!is.matrix(x)) x <- matrix(x, nrow = length(lat_idx), ncol = length(lon_idx))

  if (!is.na(missing_value)) x[x == missing_value] <- NA_real_

  lat_block <- lat[min(lat_idx):max(lat_idx)]
  lon_block <- lon[min(lon_idx):max(lon_idx)]
  lat_pos <- match(lat_keep, lat_block)
  lon_pos <- match(lon_keep, lon_block)
  x[lat_pos, lon_pos, drop = FALSE]
}

first_col <- grid_df[grid_df$c == min(grid_df$c), c("r", "lat")]
first_col <- first_col[order(first_col$r), , drop = FALSE]
lat_keep <- as.numeric(first_col$lat)

top_row <- grid_df[grid_df$r == min(grid_df$r), c("c", "lon")]
top_row <- top_row[order(top_row$c), , drop = FALSE]
lon_keep <- as.numeric(top_row$lon)

obs_sst_mat <- matrix(as.numeric(yo_full_df$Yo_full), nrow = R, ncol = C)

# Match Figure07 ensemble-panel colors.
ensemble_palette <- c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026")
yo_limits <- range(obs_sst_mat, na.rm = TRUE)

obs_site_idx <- match(as.integer(obs_df$s), as.integer(grid_df$s))
obs_points_df <- data.frame(
  X = grid_df$lon[obs_site_idx],
  Y = grid_df$lat[obs_site_idx]
)

node_df <- grid_df[, c("s", "r", "c", "lat", "lon"), drop = FALSE]
node_df$is_land <- as.logical(land_raw$is_land)
node_df <- node_df[!node_df$is_land, , drop = FALSE]
node_df$r_plot <- -node_df$r

coords <- node_df[, c("s", "r_plot", "c", "lat", "lon"), drop = FALSE]
names(coords) <- c("s", "r_plot", "c_plot", "lat", "lon")

edge_plot <- merge(edge_df, coords, by.x = "s_from", by.y = "s", all.x = TRUE, sort = FALSE)
names(edge_plot)[names(edge_plot) == "r_plot"] <- "r_from"
names(edge_plot)[names(edge_plot) == "c_plot"] <- "c_from"
names(edge_plot)[names(edge_plot) == "lat"] <- "lat_from"
names(edge_plot)[names(edge_plot) == "lon"] <- "lon_from"

edge_plot <- merge(edge_plot, coords, by.x = "s_to", by.y = "s", all.x = TRUE, sort = FALSE)
names(edge_plot)[names(edge_plot) == "r_plot"] <- "r_to"
names(edge_plot)[names(edge_plot) == "c_plot"] <- "c_to"
names(edge_plot)[names(edge_plot) == "lat"] <- "lat_to"
names(edge_plot)[names(edge_plot) == "lon"] <- "lon_to"
if (anyNA(edge_plot$r_from) || anyNA(edge_plot$r_to)) {
  stop("Some edges reference node ids not present in ocean node set.")
}

left_plot <- suppressWarnings(plot_sst_contour_map(
  Z = obs_sst_mat,
  Xp = Xp,
  Yp = Yp,
  bins = 44,
  palette = ensemble_palette,
  limits = yo_limits,
  api = NULL,
  title = "Observed SST and sampled locations",
  land_df = land_df,
  world_df = world_df
)) +
  ggplot2::geom_point(
    data = obs_points_df,
    mapping = ggplot2::aes(x = X, y = Y),
    inherit.aes = FALSE,
    color = "black",
    shape = 16,
    size = 0.42,
    alpha = 0.8
  ) +
  ggplot2::labs(fill = "SST (deg C)") +
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

right_plot <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = edge_plot,
    mapping = ggplot2::aes(x = lon_from, y = lat_from, xend = lon_to, yend = lat_to),
    inherit.aes = FALSE,
    color = "#2C7FB8",
    linewidth = 0.12,
    alpha = 0.28,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    data = node_df,
    mapping = ggplot2::aes(x = lon, y = lat),
    inherit.aes = FALSE,
    color = "#2C7FB8",
    size = 0.08,
    alpha = 0.35,
    shape = 16
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(title = "Analysis Lattice", x = NULL, y = NULL) +
  ggplot2::theme_bw(base_size = 9, base_family = "serif") +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.border = ggplot2::element_rect(color = "black", linewidth = 0.35, fill = NA),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 8.5, color = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.25, color = "black"),
    plot.title = ggplot2::element_text(hjust = 0),
    plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "pt")
  )

figure_grob <- gridExtra::arrangeGrob(
  left_plot,
  right_plot,
  nrow = 1,
  ncol = 2,
  widths = c(1.25, 1)
)

output_pdf <- file.path(fig_dir, "Figure06.pdf")
output_png <- file.path(fig_dir, "Figure06.png")

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

cropped_pdf <- tempfile(pattern = "Figure06_crop_", tmpdir = fig_dir, fileext = ".pdf")
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
  width = 10.2,
  height = 4.8,
  units = "in",
  dpi = 600,
  bg = "white"
))

message("Saved figure PDF: ", normalizePath(output_pdf, mustWork = FALSE))
message("Saved figure PNG: ", normalizePath(output_png, mustWork = FALSE))
