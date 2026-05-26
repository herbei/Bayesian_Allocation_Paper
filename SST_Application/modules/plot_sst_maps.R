# SST plotting helpers: UQ26 contour style + geographic boundaries.

source_pacific_plot_utils <- function(plot_utils_path = file.path("..", "R", "modules", "plot_utils.R")) {
  if (!file.exists(plot_utils_path)) {
    stop("Could not find UQ26 plot utility file: ", plot_utils_path)
  }
  source(plot_utils_path, local = parent.frame())
}

to_0360_lon <- function(lon) {
  ifelse(lon < 0, lon + 360, lon)
}

build_world_boundaries <- function(lon_limits, lat_limits, lon_buffer = 2, lat_buffer = 2,
                                   jump_threshold = 30) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for map boundaries.")
  }
  if (!requireNamespace("maps", quietly = TRUE)) {
    stop("Package 'maps' is required for map boundaries.")
  }

  world <- ggplot2::map_data("world")
  world$long <- to_0360_lon(world$long)

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
    stop("Package 'maps' is required for land masking.")
  }

  x <- as.vector(Xp)
  y <- as.vector(Yp)
  x180 <- ifelse(x > 180, x - 360, x)
  region <- maps::map.where("world", x = x180, y = y)

  data.frame(
    X = x,
    Y = y,
    is_land = !is.na(region)
  )
}

default_non_pacific_rectangles <- function() {
  # Connected boxes in 0..360 longitude coordinates.
  # These target Gulf of Mexico, Caribbean, and western Atlantic only.
  data.frame(
    region = c(
      "gulf_mexico",
      "gulf_carib_connector",
      "caribbean",
      "atlantic_tropical",
      "atlantic_north"
    ),
    lon_min = c(258.5, 258.5, 264.5, 280.5, 276.5),
    lon_max = c(281.5, 263.5, 286.5, 288.5, 288.5),
    lat_min = c(18.5,  15.5,   9.5,  -4.5,  23.5),
    lat_max = c(31.5,  17.5,  23.5,  23.5,  44.5),
    stringsAsFactors = FALSE
  )
}

build_pacific_ocean_mask <- function(Xp, Yp, land_df = NULL, rectangles = NULL) {
  if (is.null(land_df)) {
    land_df <- build_land_mask(Xp = Xp, Yp = Yp)
  }
  if (is.null(rectangles)) {
    rectangles <- default_non_pacific_rectangles()
  }

  x_vec <- as.vector(Xp)
  y_vec <- as.vector(Yp)
  ocean_vec <- !land_df$is_land
  rect_mask <- rep(FALSE, length(x_vec))
  for (i in seq_len(nrow(rectangles))) {
    rect_mask <- rect_mask | (
      x_vec >= rectangles$lon_min[i] & x_vec <= rectangles$lon_max[i] &
        y_vec >= rectangles$lat_min[i] & y_vec <= rectangles$lat_max[i]
    )
  }

  data.frame(
    X = x_vec,
    Y = y_vec,
    is_non_pacific_ocean = ocean_vec & rect_mask
  )
}

add_non_pacific_mask <- function(p, ocean_mask_df, fill = "grey80") {
  if (is.null(ocean_mask_df) || nrow(ocean_mask_df) == 0) return(p)

  ocean_cells <- ocean_mask_df[ocean_mask_df$is_non_pacific_ocean, , drop = FALSE]
  if (nrow(ocean_cells) == 0) return(p)

  p + ggplot2::geom_tile(
    data = ocean_cells,
    mapping = ggplot2::aes(x = X, y = Y),
    inherit.aes = FALSE,
    width = 1,
    height = 1,
    fill = fill,
    alpha = 1
  )
}

add_land_mask <- function(p, land_df, fill = "grey80") {
  if (is.null(land_df) || nrow(land_df) == 0) return(p)

  land_cells <- land_df[land_df$is_land, , drop = FALSE]
  if (nrow(land_cells) == 0) return(p)

  p + ggplot2::geom_tile(
    data = land_cells,
    mapping = ggplot2::aes(x = X, y = Y),
    inherit.aes = FALSE,
    width = 1,
    height = 1,
    fill = fill,
    alpha = 1
  )
}

add_world_boundaries <- function(p, world_df, color = "black", linewidth = 0.25, alpha = 0.85) {
  if (is.null(world_df) || nrow(world_df) == 0) return(p)

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

plot_sst_contour_map <- function(Z, Xp, Yp, bins = 40, palette, api = NULL, zero_white = FALSE,
                                 limits = NULL, title = NULL, land_df = NULL, world_df = NULL,
                                 ocean_mask_df = NULL, land_fill = "grey80", ocean_fill = "grey80",
                                 obs_point_size = 0.35, obs_point_alpha = 0.75,
                                 obs_point_shape = 46) {
  p <- plot_filled_contour(
    Z = Z,
    Xp = Xp,
    Yp = Yp,
    bins = bins,
    palette = palette,
    api = NULL,
    zero_white = zero_white,
    limits = limits
  )

  if (!is.null(title)) {
    p <- p + ggplot2::labs(title = title)
  }

  p <- add_non_pacific_mask(p, ocean_mask_df = ocean_mask_df, fill = ocean_fill)
  p <- add_land_mask(p, land_df = land_df, fill = land_fill)
  p <- add_world_boundaries(p, world_df = world_df)

  if (!is.null(api)) {
    api <- as.integer(api)
    grid_df <- data.frame(
      idx = seq_along(as.vector(Xp)),
      X = as.vector(Xp),
      Y = as.vector(Yp)
    )
    api <- api[is.finite(api) & api >= 1L & api <= nrow(grid_df)]
    if (length(api) > 0) {
      p <- p + ggplot2::geom_point(
        data = grid_df[api, , drop = FALSE],
        mapping = ggplot2::aes(x = X, y = Y),
        inherit.aes = FALSE,
        shape = obs_point_shape,
        size = obs_point_size,
        color = "black",
        alpha = obs_point_alpha
      )
    }
  }

  p
}
