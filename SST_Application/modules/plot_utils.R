plot_filled_contour <- function(Z, Xp, Yp, bins = 40, palette, api = NULL, zero_white = FALSE, limits = NULL) {
  plot_title <- deparse(substitute(Z))

  df <- data.frame(
    X = as.vector(Xp),
    Y = as.vector(Yp),
    Z = as.vector(Z)
  )

  grid_df <- data.frame(
    idx = seq_along(as.vector(Xp)),
    X = as.vector(Xp),
    Y = as.vector(Yp)
  )

  fill_values <- NULL
  z_min <- if (is.null(limits)) min(df$Z, na.rm = TRUE) else limits[1]
  z_max <- if (is.null(limits)) max(df$Z, na.rm = TRUE) else limits[2]

  if (zero_white && z_min < 0 && z_max > 0) {
    n_pal <- length(palette)
    mid_idx <- ceiling(n_pal / 2)
    zero_pos <- (0 - z_min) / (z_max - z_min)
    lower_vals <- seq(0, zero_pos, length.out = mid_idx)
    upper_vals <- seq(zero_pos, 1, length.out = n_pal - mid_idx + 1)
    fill_values <- c(lower_vals, upper_vals[-1])
  }

  p <- ggplot(df, aes(x = X, y = Y, z = Z)) +
    geom_contour_filled(
      aes(fill = after_stat((level_low + level_high) / 2)),
      bins = bins
    ) +
    coord_equal(expand = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colors = palette, values = fill_values, limits = limits, oob = scales::squish) +
    labs(title = plot_title, x = NULL, y = NULL, fill = NULL) +
    theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(0, 0, 2, 0, unit = "pt"))
    )

  if (is.null(api) && exists("api", envir = parent.frame(), inherits = TRUE)) {
    api <- get("api", envir = parent.frame(), inherits = TRUE)
  }

  if (!is.null(api)) {
    api <- as.integer(api)
    api <- api[is.finite(api) & api >= 1L & api <= nrow(grid_df)]

    p <- p + geom_point(
      data = grid_df[api, , drop = FALSE],
      aes(x = X, y = Y),
      inherit.aes = FALSE,
      size = 0.8,
      color = "black",
      alpha = 0.75
    )
  }

  p
}

redblue <- function(m) {
  if (m %% 2 == 0) {
    m1 <- m * 0.5
    r <- (0:(m1 - 1)) / max(m1 - 1, 1)
    g <- r
    r <- c(r, rep(1, m1))
    g <- c(g, rev(g))
    b <- rev(r)
  } else {
    m1 <- floor(m * 0.5)
    r <- (0:(m1 - 1)) / max(m1, 1)
    g <- r
    r <- c(r, rep(1, m1 + 1))
    g <- c(g, 1, rev(g))
    b <- rev(r)
  }
  cbind(r, g, b)
}

build_plot_palettes <- function() {
  rb <- redblue(256)
  red_blue <- rgb(rb[, 1], rb[, 2], rb[, 3], maxColorValue = 1)
  my_palette <- rev(c(
    "#5e1812", "#a02410", "#ce3118", "#e33f1c", "#e66222",
    "#e98e2a", "#eebc35", "#f5ed40", "#d4f64e", "#abf46f",
    "#89f299", "#6ff1c5", "#63ecf4", "#4dbcf1", "#388bee",
    "#2358ef", "#0e1eec", "#0400de", "#0300ad", "#07025b"
  ))
  list(rb = rb, red_blue = red_blue, my_palette = my_palette)
}
