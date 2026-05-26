#!/usr/bin/env Rscript

# Figure 07: publication-ready 2x2 panel of four E3SM ensemble members.

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
fig_dir <- resolve_figure_output_dir(figure_root, run_label = "static")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ensemble_rds <- file.path("Data", "processed", "pacific_sst_ensemble_subindexed_201201.rds")
if (!file.exists(ensemble_rds)) {
  stop("Missing ensemble data file: ", ensemble_rds)
}

ens <- readRDS(ensemble_rds)
required_fields <- c("meta", "domain", "members", "fields")
if (!all(required_fields %in% names(ens))) {
  stop(
    "Ensemble RDS is missing required fields. Expected: ",
    paste(required_fields, collapse = ", ")
  )
}

target_members <- c("r01", "r07", "r14", "r21")
member_idx <- match(target_members, ens$members)
if (anyNA(member_idx)) {
  missing_members <- target_members[is.na(member_idx)]
  stop("Requested members not found in ensemble object: ", paste(missing_members, collapse = ", "))
}

lat_keep <- as.numeric(ens$domain$lat)
lon_keep <- as.numeric(ens$domain$lon)
R <- length(lat_keep)
C <- length(lon_keep)

Xp <- matrix(rep(lon_keep, each = R), nrow = R, ncol = C)
Yp <- matrix(rep(lat_keep, times = C), nrow = R, ncol = C)

land_file <- file.path("Data", "processed", "pacific_land_mask.csv")
if (file.exists(land_file)) {
  land_raw <- utils::read.csv(land_file, stringsAsFactors = FALSE)
  land_df <- data.frame(
    X = as.numeric(land_raw$lon),
    Y = as.numeric(land_raw$lat),
    is_land = as.logical(land_raw$is_land)
  )
} else {
  land_df <- build_land_mask(Xp = Xp, Yp = Yp)
}

world_df <- build_world_boundaries(
  lon_limits = range(Xp, na.rm = TRUE),
  lat_limits = range(Yp, na.rm = TRUE)
)

palette_nondiv <- c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026")
common_limits <- range(ens$fields, na.rm = TRUE)

extract_legend <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(idx) == 0L) return(NULL)
  g$grobs[[idx[1]]]
}

panel_theme <- ggplot2::theme_bw(base_size = 9, base_family = "serif") +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = "black", linewidth = 0.35),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 8, color = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.25),
    plot.title = ggplot2::element_text(hjust = 0),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8, color = "black"),
    legend.key.height = grid::unit(7, "pt"),
    legend.key.width = grid::unit(30, "pt"),
    plot.margin = ggplot2::margin(1, 3, 1, 3, unit = "pt")
  )

panel_plots <- vector("list", length(target_members))
for (i in seq_along(target_members)) {
  panel <- suppressWarnings(plot_sst_contour_map(
    Z = ens$fields[, , member_idx[i]],
    Xp = Xp,
    Yp = Yp,
    bins = 45,
    palette = palette_nondiv,
    limits = common_limits,
    title = sprintf("(%s) Member %s", LETTERS[i], target_members[i]),
    land_df = land_df,
    world_df = world_df
  )) +
    ggplot2::labs(fill = "SST (deg C)") +
    panel_theme

  if (i <= 2L) {
    panel <- panel + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(1, 3, 0, 3, unit = "pt")
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

legend_plot <- panel_plots[[1]] +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = grid::unit(3.4, "in"),
      barheight = grid::unit(0.14, "in"),
      ticks.colour = "black",
      frame.colour = "black",
      title.theme = ggplot2::element_text(size = 9, color = "black"),
      label.theme = ggplot2::element_text(size = 8, color = "black")
    )
  )
legend_grob <- suppressWarnings(extract_legend(legend_plot))
if (is.null(legend_grob)) {
  stop("Could not extract legend from panel plot.")
}
legend_height <- sum(legend_grob$heights) + grid::unit(4, "pt")

panel_plots_nolegend <- lapply(
  panel_plots,
  function(p) p + ggplot2::theme(legend.position = "none")
)

panel_grobs <- lapply(panel_plots_nolegend, ggplot2::ggplotGrob)
common_widths <- do.call(grid::unit.pmax, lapply(panel_grobs, function(g) g$widths))
panel_grobs <- lapply(panel_grobs, function(g) {
  g$widths <- common_widths
  g
})

top_row_grob <- suppressWarnings(gridExtra::arrangeGrob(
  grobs = panel_grobs[1:2],
  ncol = 2,
  padding = grid::unit(0, "pt")
))

bottom_row_grob <- suppressWarnings(gridExtra::arrangeGrob(
  grobs = panel_grobs[3:4],
  ncol = 2,
  padding = grid::unit(0, "pt")
))

figure_grob <- gridExtra::arrangeGrob(
  top_row_grob,
  bottom_row_grob,
  legend_grob,
  ncol = 1,
  heights = grid::unit.c(grid::unit(1, "null"), grid::unit(1, "null"), legend_height)
)

output_pdf <- file.path(fig_dir, "Figure07.pdf")
output_png <- file.path(fig_dir, "Figure07.png")

pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
suppressWarnings(ggplot2::ggsave(
  filename = output_pdf,
  plot = figure_grob,
  device = pdf_device,
  width = 7.5,
  height = 5.3,
  units = "in",
  bg = "white"
))

cropped_pdf <- tempfile(pattern = "Figure07_crop_", tmpdir = fig_dir, fileext = ".pdf")
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
  width = 7.5,
  height = 5.3,
  units = "in",
  dpi = 600,
  bg = "white"
))

message("Saved figure PDF: ", normalizePath(output_pdf, mustWork = FALSE))
message("Saved figure PNG: ", normalizePath(output_png, mustWork = FALSE))
