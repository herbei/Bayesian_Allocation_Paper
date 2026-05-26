#!/usr/bin/env Rscript

# Figure 09b: posterior histograms for sigma^2_yo and sigma^2_yc
# from the selected SST ensemble run.

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
    if (file.exists(file.path(cand, "Data", "results"))) {
      return(cand)
    }
  }
  stop("Could not locate project root containing Data/results.")
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

source(file.path(project_root, "Figures", "src", "figure_io_utils.R"))

figure_root <- file.path(project_root, "Figures")
best_rdata <- resolve_figure_chain_file(figure_root, "SST_FIGURE09B_CHAIN_FILE")
fig_dir <- resolve_figure_output_dir(figure_root, chain_file = best_rdata)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(best_rdata)) {
  stop("Missing best-run RData file: ", best_rdata)
}

result_env <- new.env(parent = emptyenv())
load(best_rdata, envir = result_env)
if (!exists("mcmc_artifact", envir = result_env, inherits = FALSE)) {
  stop("Best-run file does not contain mcmc_artifact: ", best_rdata)
}

artifact <- result_env$mcmc_artifact
store <- artifact$output$store
run_config <- artifact$run_config
variance_update_start <- if (!is.null(run_config$update_sigma2_after)) {
  as.integer(run_config$update_sigma2_after)
} else {
  15000L
}
iter_chain <- as.integer(store$iter_chain)
keep_idx <- which(iter_chain > variance_update_start)

sig2o_post <- as.numeric(store$sig2o_chain[keep_idx])
sig2c_post <- as.numeric(store$sig2c_chain[keep_idx])

if (length(sig2o_post) == 0L || length(sig2c_post) == 0L) {
  stop("No sigma draws found after the variance-update threshold in best-run artifact.")
}

hist_theme <- ggplot2::theme_bw(base_size = 9, base_family = "serif") +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = "black", linewidth = 0.35),
    axis.text = ggplot2::element_text(size = 8, color = "black"),
    axis.title = ggplot2::element_text(size = 9),
    axis.ticks = ggplot2::element_line(linewidth = 0.25),
    plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0),
    plot.margin = ggplot2::margin(3, 4, 3, 4, unit = "pt")
  )

sig2o_df <- data.frame(value = sig2o_post)
sig2c_df <- data.frame(value = sig2c_post)
sig2o_mean <- mean(sig2o_post)
sig2c_mean <- mean(sig2c_post)
sig2o_label <- sprintf("mean = %.3f", sig2o_mean)
sig2c_label <- sprintf("mean = %.3f", sig2c_mean)

p_sig2o <- ggplot(sig2o_df, aes(x = value)) +
  geom_histogram(
    bins = 20,
    fill = "#1b9e77",
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = sig2o_mean,
    color = "black",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  annotate(
    geom = "label",
    x = Inf,
    y = Inf,
    hjust = 1.04,
    vjust = 1.10,
    label = sig2o_label,
    size = 2.9,
    family = "mono",
    label.size = 0.2,
    fill = grDevices::adjustcolor("white", alpha.f = 0.9)
  ) +
  labs(
    title = expression("(A) Posterior of " * sigma[yo]^2),
    x = expression(sigma[yo]^2),
    y = "Count"
  ) +
  hist_theme

p_sig2c <- ggplot(sig2c_df, aes(x = value)) +
  geom_histogram(
    bins = 20,
    fill = "#d95f02",
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = sig2c_mean,
    color = "black",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  annotate(
    geom = "label",
    x = Inf,
    y = Inf,
    hjust = 1.04,
    vjust = 1.10,
    label = sig2c_label,
    size = 2.9,
    family = "mono",
    label.size = 0.2,
    fill = grDevices::adjustcolor("white", alpha.f = 0.9)
  ) +
  labs(
    title = expression("(B) Posterior of " * sigma[yc]^2),
    x = expression(sigma[yc]^2),
    y = "Count"
  ) +
  hist_theme

figure_grob <- gridExtra::arrangeGrob(
  p_sig2o,
  p_sig2c,
  nrow = 1,
  ncol = 2
)

output_pdf <- file.path(fig_dir, "Figure09b.pdf")
output_png <- file.path(fig_dir, "Figure09b.png")

figure_width_in <- 7.4
figure_height_in <- 2.64
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

cropped_pdf <- tempfile(pattern = "Figure09b_crop_", tmpdir = fig_dir, fileext = ".pdf")
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
