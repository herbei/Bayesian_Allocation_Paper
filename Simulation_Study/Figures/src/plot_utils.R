build_world_boundaries <- function(lon_limits, lat_limits, lon_buffer = 2, lat_buffer = 2,
                                   jump_threshold = 30) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for contour plotting.", call. = FALSE)
  }
  if (!requireNamespace("maps", quietly = TRUE)) {
    stop("Package `maps` is required for coastline overlays.", call. = FALSE)
  }

  world <- ggplot2::map_data("world")
  world$seg_group <- NA_character_
  groups <- unique(world$group)
  for (g in groups) {
    idx <- which(world$group == g)
    lon_g <- world$long[idx]
    split_id <- cumsum(c(0L, as.integer(abs(diff(lon_g)) > jump_threshold)))
    world$seg_group[idx] <- paste0(g, "_", split_id)
  }

  lon_min <- min(lon_limits) - lon_buffer
  lon_max <- max(lon_limits) + lon_buffer
  lat_min <- min(lat_limits) - lat_buffer
  lat_max <- max(lat_limits) + lat_buffer

  keep <- world$long >= lon_min & world$long <= lon_max &
    world$lat >= lat_min & world$lat <= lat_max

  world[keep, c("long", "lat", "seg_group")]
}

build_land_mask <- function(Xp, Yp) {
  if (!requireNamespace("maps", quietly = TRUE)) {
    stop("Package `maps` is required for land masking.", call. = FALSE)
  }

  x <- as.vector(Xp)
  y <- as.vector(Yp)
  region <- maps::map.where("world", x = x, y = y)

  data.frame(
    X = x,
    Y = y,
    is_land = !is.na(region)
  )
}

add_land_mask <- function(p, land_df, width, height, fill = "grey80") {
  if (is.null(land_df) || nrow(land_df) == 0L) {
    return(p)
  }

  land_cells <- land_df[land_df$is_land, , drop = FALSE]
  if (nrow(land_cells) == 0L) {
    return(p)
  }

  p + ggplot2::geom_tile(
    data = land_cells,
    mapping = ggplot2::aes(x = X, y = Y),
    inherit.aes = FALSE,
    width = width,
    height = height,
    fill = fill,
    alpha = 1
  )
}

add_world_boundaries <- function(p, world_df, color = "black", linewidth = 0.25, alpha = 0.85) {
  if (is.null(world_df) || nrow(world_df) == 0L) {
    return(p)
  }

  p + ggplot2::geom_path(
    data = world_df,
    mapping = ggplot2::aes(x = long, y = lat, group = seg_group),
    inherit.aes = FALSE,
    color = color,
    linewidth = linewidth,
    alpha = alpha,
    lineend = "round"
  )
}

build_oxy_palette <- function() {
  c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026")
}

build_oxy_panel_theme <- function() {
  ggplot2::theme_bw(base_size = 9, base_family = "serif") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.35),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 8, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25),
      plot.title = ggplot2::element_text(
        size = 10,
        face = "bold",
        hjust = 0,
        margin = ggplot2::margin(0, 0, 1, 0, unit = "pt")
      ),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8, color = "black"),
      legend.key.height = grid::unit(7, "pt"),
      legend.key.width = grid::unit(30, "pt"),
      plot.margin = ggplot2::margin(1, 3, 1, 3, unit = "pt")
    )
}

plot_oxy_contour <- function(Z, Xp, Yp, title, fill_label = "OXY", land_df = NULL, world_df = NULL,
                             x_limits = NULL, y_limits = NULL, fill_limits = NULL, bins = 45L) {
  df <- data.frame(
    X = as.vector(Xp),
    Y = as.vector(Yp),
    Z = as.vector(Z)
  )
  x_step <- if (length(unique(df$X)) > 1L) median(diff(sort(unique(df$X)))) else 1
  y_step <- if (length(unique(df$Y)) > 1L) median(diff(sort(unique(df$Y)))) else 1

  if (is.null(fill_limits)) {
    fill_limits <- range(df$Z, na.rm = TRUE)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = X, y = Y, z = Z)) +
    ggplot2::geom_contour_filled(
      ggplot2::aes(fill = after_stat((level_low + level_high) / 2)),
      bins = bins
    ) +
    ggplot2::coord_equal(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = x_limits) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = y_limits) +
    ggplot2::scale_fill_gradientn(
      colors = build_oxy_palette(),
      limits = fill_limits,
      oob = scales::squish
    ) +
    ggplot2::labs(title = title, fill = fill_label) +
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
    ) +
    build_oxy_panel_theme()

  p <- add_land_mask(p, land_df = land_df, width = x_step, height = y_step)
  add_world_boundaries(p, world_df = world_df)
}

extract_legend <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1L)) == "guide-box")
  if (length(idx) == 0L) {
    return(NULL)
  }
  g$grobs[[idx[1L]]]
}

crop_pdf_in_place <- function(pdf_path) {
  if (!file.exists(pdf_path)) {
    stop("Cannot crop missing PDF: ", pdf_path, call. = FALSE)
  }

  pdfcrop_bin <- Sys.which("pdfcrop")
  if (!nzchar(pdfcrop_bin)) {
    stop("`pdfcrop` is not available in PATH.", call. = FALSE)
  }

  cropped_pdf <- tempfile(
    pattern = paste0(tools::file_path_sans_ext(basename(pdf_path)), "_crop_"),
    tmpdir = dirname(pdf_path),
    fileext = ".pdf"
  )

  crop_status <- suppressWarnings(system2(
    pdfcrop_bin,
    args = c(pdf_path, cropped_pdf),
    stdout = FALSE,
    stderr = FALSE
  ))

  if (!identical(crop_status, 0L) || !file.exists(cropped_pdf)) {
    if (file.exists(cropped_pdf)) {
      unlink(cropped_pdf)
    }
    stop("`pdfcrop` failed for ", pdf_path, call. = FALSE)
  }

  ok_rename <- file.rename(cropped_pdf, pdf_path)
  if (!isTRUE(ok_rename)) {
    ok_remove <- file.remove(pdf_path)
    if (!isTRUE(ok_remove)) {
      unlink(cropped_pdf)
      stop("`pdfcrop` succeeded but failed to replace the original PDF.", call. = FALSE)
    }
    ok_rename <- file.rename(cropped_pdf, pdf_path)
  }

  if (!isTRUE(ok_rename)) {
    if (file.exists(cropped_pdf)) {
      unlink(cropped_pdf)
    }
    stop("`pdfcrop` succeeded but failed to replace the original PDF.", call. = FALSE)
  }

  invisible(normalizePath(pdf_path, mustWork = FALSE))
}

save_panel_figure <- function(panel_plots, output_pdf, output_png) {
  if (length(panel_plots) != 4L) {
    stop("`panel_plots` must contain exactly 4 panels.", call. = FALSE)
  }

  legend_plot <- panel_plots[[1L]]
  legend_grob <- suppressWarnings(extract_legend(legend_plot))
  if (is.null(legend_grob)) {
    stop("Could not extract legend from panel plot.", call. = FALSE)
  }
  legend_height <- sum(legend_grob$heights) + grid::unit(4, "pt")

  panel_plots_nolegend <- lapply(
    panel_plots,
    function(p) p + ggplot2::theme(legend.position = "none")
  )

  panel_grobs <- lapply(panel_plots_nolegend, ggplot2::ggplotGrob)
  common_widths <- do.call(grid::unit.pmax, lapply(panel_grobs, function(g) g$widths))
  common_heights <- do.call(grid::unit.pmax, lapply(panel_grobs, function(g) g$heights))
  panel_grobs <- lapply(panel_grobs, function(g) {
    g$widths <- common_widths
    g$heights <- common_heights
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
  crop_pdf_in_place(output_pdf)

  suppressWarnings(ggplot2::ggsave(
    filename = output_png,
    plot = figure_grob,
    width = 7.5,
    height = 5.3,
    units = "in",
    dpi = 600,
    bg = "white"
  ))
}

resolve_figure_chain_file <- function(fig_root,
                                      env_names = character(0),
                                      default_filename = file.path("..", "Data", "results", "simulation_study_selected_run.RData")) {
  env_names <- unique(c(env_names, "OXY_FIGURE_CHAIN_FILE"))
  env_names <- env_names[nzchar(env_names)]

  for (env_name in env_names) {
    env_value <- Sys.getenv(env_name, unset = "")
    if (nzchar(env_value)) {
      return(normalizePath(env_value, winslash = "/", mustWork = FALSE))
    }
  }

  normalizePath(file.path(fig_root, default_filename), winslash = "/", mustWork = FALSE)
}

figure_run_label <- function(chain_file) {
  label <- basename(chain_file)
  label <- sub("(?i)\\.rdata$", "", label, perl = TRUE)
  label <- gsub("[^A-Za-z0-9._-]+", "_", label)
  if (!nzchar(label)) {
    label <- "run"
  }
  label
}

resolve_figure_output_dir <- function(fig_root, chain_file = NULL) {
  output_dir_env <- Sys.getenv("OXY_FIGURE_OUTPUT_DIR", unset = "")
  if (nzchar(output_dir_env)) {
    return(normalizePath(output_dir_env, winslash = "/", mustWork = FALSE))
  }

  if (!is.null(chain_file) && nzchar(chain_file)) {
    return(normalizePath(
      file.path(fig_root, "out", figure_run_label(chain_file)),
      winslash = "/",
      mustWork = FALSE
    ))
  }

  normalizePath(fig_root, winslash = "/", mustWork = FALSE)
}
