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
  library(gridExtra)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

source(file.path(project_root, "Figures", "src", "plot_utils.R"))

figure_root <- file.path(project_root, "Figures")
chain_file <- resolve_figure_chain_file(
  figure_root,
  "OXY_FIGURE05_CHAIN_FILE"
)
fig_dir <- resolve_figure_output_dir(figure_root, chain_file)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(chain_file)) {
  stop("Missing chain file: ", chain_file, call. = FALSE)
}

result_env <- new.env(parent = emptyenv())
load(chain_file, envir = result_env)
if (!exists("mcmc_artifact", envir = result_env, inherits = FALSE)) {
  stop("Best-run file does not contain mcmc_artifact: ", chain_file, call. = FALSE)
}

artifact <- result_env$mcmc_artifact
store <- artifact$output$store
keep_idx <- as.integer(artifact$output$keep_idx %||% integer(0))
B_chain <- store$B_chain %||% NULL

if (is.null(B_chain) || length(dim(B_chain)) != 2L || ncol(B_chain) == 0L) {
  stop("No B-chain columns available for PIP diagnostics in selected artifact.", call. = FALSE)
}
if (length(keep_idx) == 0L) {
  stop("No kept post-burn draws available in selected artifact.", call. = FALSE)
}

B_post <- B_chain[keep_idx, , drop = FALSE]
k_vals <- as.integer(sub("^k", "", colnames(B_post)))
pip_X <- colMeans(B_post)
pip_beta <- 1 - pip_X
pip_plot_df <- rbind(
  data.frame(k = k_vals, component = "X", prob = pip_X),
  data.frame(k = k_vals, component = "beta", prob = pip_beta)
)
pip_plot_df$component <- factor(pip_plot_df$component, levels = c("beta", "X"))

display_threshold <- 0.5
focus_threshold <- 0.25
focus_cutoff <- display_threshold + focus_threshold

# Keep the grid structure visible while reserving the color scale for shown cells.
pip_full_df <- pip_plot_df
pip_full_df$prob_display <- ifelse(pip_full_df$prob > display_threshold, pip_full_df$prob, NA_real_)

pip_focus_df <- pip_plot_df
pip_focus_df$prob_display <- ifelse(pip_focus_df$prob > focus_cutoff, pip_focus_df$prob, NA_real_)

pip_theme <- ggplot2::theme_bw(base_size = 9, base_family = "serif") +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = "black", linewidth = 0.35),
    axis.text = ggplot2::element_text(size = 8, color = "black"),
    axis.title = ggplot2::element_text(size = 9),
    axis.ticks = ggplot2::element_line(linewidth = 0.25),
    plot.margin = ggplot2::margin(3, 4, 3, 4, unit = "pt")
  )

extract_legend <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1L)) == "guide-box")
  if (length(idx) == 0L) return(NULL)
  g$grobs[[idx[1L]]]
}

compute_x_breaks <- function(k_vals, step = 50L) {
  k_vals <- sort(unique(as.integer(k_vals)))
  lower <- ceiling(min(k_vals) / step) * step
  upper <- floor(max(k_vals) / step) * step

  if (lower <= upper) {
    seq(lower, upper, by = step)
  } else {
    unique(as.integer(round(seq(min(k_vals), max(k_vals), length.out = 5L))))
  }
}

bin_breaks <- seq(display_threshold, 1, by = 0.1)
bin_palette <- c(
  "#fee8c8",
  "#fdbb84",
  "#fc8d59",
  "#e34a33",
  "#b30000"
)

build_pip_heatmap <- function(df, title_text = NULL, show_legend = TRUE, breaks = NULL, labels = NULL) {
  x_scale <- ggplot2::scale_x_continuous(
    expand = c(0, 0),
    breaks = breaks %||% compute_x_breaks(df$k)
  )
  if (!is.null(labels)) {
    x_scale <- ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = breaks %||% compute_x_breaks(df$k),
      labels = labels
    )
  }

  p <- ggplot(df, aes(x = k, y = component, fill = prob_display)) +
    geom_tile(width = 1, height = 0.82) +
    x_scale +
    scale_y_discrete(
      expand = c(0, 0),
      labels = c(expression(beta), expression(X))
    ) +
    scale_fill_stepsn(
      colours = bin_palette,
      limits = c(display_threshold, 1),
      breaks = bin_breaks,
      labels = formatC(bin_breaks, format = "f", digits = 1),
      na.value = "#ececec",
      name = "Probability"
    ) +
    guides(
      fill = ggplot2::guide_coloursteps(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(2.4, "in"),
        barheight = grid::unit(0.14, "in")
      )
    ) +
    labs(
      title = title_text,
      x = "Mode index k",
      y = NULL
    ) +
    pip_theme
  if (show_legend) {
    p <- p + ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(size = 8.5),
      legend.text = ggplot2::element_text(size = 7.5),
      legend.key.width = grid::unit(42, "pt"),
      legend.key.height = grid::unit(8, "pt")
    )
  } else {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p
}

full_breaks <- compute_x_breaks(k_vals)

probe_pdf <- tempfile(fileext = ".pdf")
grDevices::pdf(file = probe_pdf, width = 7.4, height = 4.1)
legend_plot <- build_pip_heatmap(
  df = pip_full_df,
  title_text = NULL,
  show_legend = TRUE,
  breaks = full_breaks
)

legend_grob <- suppressWarnings(extract_legend(legend_plot))
if (is.null(legend_grob)) {
  stop("Could not extract legend for Figure05.")
}

p_pip_full <- build_pip_heatmap(
  df = pip_full_df,
  title_text = "(A) Posterior allocation probabilities > 0.5 across all modes",
  show_legend = FALSE,
  breaks = full_breaks
)

p_pip_focus <- build_pip_heatmap(
  df = pip_focus_df,
  title_text = sprintf("(B) Highlighted modes with posterior allocation probability > %.2f", focus_cutoff),
  show_legend = FALSE,
  breaks = full_breaks
)

panel_grob <- gridExtra::arrangeGrob(
  p_pip_full,
  p_pip_focus,
  ncol = 1
)

figure_grob <- gridExtra::arrangeGrob(
  legend_grob,
  panel_grob,
  ncol = 1,
  heights = c(0.12, 1)
)
grDevices::dev.off()
unlink(probe_pdf)

output_pdf <- file.path(fig_dir, "Figure05.pdf")
output_png <- file.path(fig_dir, "Figure05.png")

figure_width_in <- 7.4
figure_height_in <- 4.1
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

message("Chain file: ", normalizePath(chain_file, mustWork = FALSE))
message("Saved figure PDF: ", normalizePath(output_pdf, mustWork = FALSE))
message("Saved figure PNG: ", normalizePath(output_png, mustWork = FALSE))
