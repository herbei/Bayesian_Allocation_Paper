# UQ26: Sea Surface Temperature (SST) real-data pipeline
#
# This script downloads observation-based gridded SST (OISST v2 monthly mean)
# from NOAA PSL THREDDS and builds observation/grid inputs compatible with the
# manuscript's lattice setup:
#   - Regular R x C rectangular lattice
#   - Physical observations Yo at a subset of sites O
#   - Column-major site indexing: s = r + (c-1)R
#
# E3SM ensemble downloads are handled by download_e3sm_ensembles.R, and E3SM
# extraction is handled by data_retrieval_ens.R.
#
# Required packages: ncdf4, Matrix, maps

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
sst_app_dir <- if (basename(script_dir) == "DO_THIS_FIRST") {
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
} else {
  script_dir
}
setwd(sst_app_dir)

project_lib <- file.path(sst_app_dir, ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(ncdf4)
  library(Matrix)
})

# ---------------------------
# 0) User-configurable inputs
# ---------------------------

# Output folder (relative to the SST_Application directory)
out_dir <- "Data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dir.create(file.path(out_dir, "raw"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "processed"), showWarnings = FALSE, recursive = TRUE)

# Month to extract (YYYY-MM). Keep this aligned with data_retrieval_ens.R.
target_ym <- "2012-01"

# Lattice dimensions (legacy manuscript default)
# R <- 19
# C <- 37
# S <- R * C

# Domain limits (legacy tropical Pacific window)
# lat_keep <- seq(  8.5,  -9.5, by = -1)  # 19 values, north-to-south
# lon_keep <- seq(150.5, 186.5, by =  1)  # 37 values, west-to-east

# Pacific-wide window on the 0..360 longitude grid.
# OISST/E3SM regridded files use 1-degree centers:
#   lon = 0.5, 1.5, ..., 359.5; lat = 89.5, 88.5, ..., -89.5.
# Requested limits: lat -45..45, lon 130..289.
# On this grid, nearest centers are lat -44.5..44.5 and lon 130.5..288.5.
lat_keep <- seq( 44.5, -44.5, by = -1)
lon_keep <- seq(130.5, 288.5, by =  1)

R <- length(lat_keep)
C <- length(lon_keep)
S <- R * C

# Number of observed sites (subset size). Berliner et al. synthetic example used 335.
N_obs <- 1000
set.seed(1)  # reproducible observed-site selection

# ---------------------------
# 1) Download raw observation netCDF file
# ---------------------------

# (a) Observations: NOAA OISST v2 monthly mean (1x1)
# NOAA PSL THREDDS HTTPServer file download
url_obs <- "https://psl.noaa.gov/thredds/fileServer/Datasets/noaa.oisst.v2/sst.mnmean.nc"
file_obs <- file.path(out_dir, "raw", "oisst_v2_sst_mnmean.nc")

# Helper: download if missing
maybe_download <- function(url, dest) {
  if (!file.exists(dest)) {
    message("Downloading: ", url)
    utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  } else {
    message("Already present: ", dest)
  }
}

maybe_download(url_obs, file_obs)

# ---------------------------
# 2) Utilities
# ---------------------------

# Convert OISST time (days since 1800-01-01) to year-month strings "YYYY-MM"
# OISST is monthly means. We use the mid-month date for matching.
get_ym_oisst <- function(nc) {
  time <- ncvar_get(nc, "time")
  origin <- as.Date("1800-01-01")
  dates <- origin + time
  format(dates, "%Y-%m")
}

# Extract a (lat,lon) window at a single time index and return as a matrix [R x C]
# with rows corresponding to lat_keep order and columns corresponding to lon_keep order.
extract_window <- function(nc, varname, t_index, lat_keep, lon_keep, missing_value = NA_real_,
                           apply_scale = TRUE) {
  lat <- ncvar_get(nc, "lat")
  lon <- ncvar_get(nc, "lon")

  lat_idx <- match(lat_keep, lat)
  lon_idx <- match(lon_keep, lon)

  if (any(is.na(lat_idx))) stop("Some lat_keep values not found in dataset lat array.")
  if (any(is.na(lon_idx))) stop("Some lon_keep values not found in dataset lon array.")

  v <- nc$var[[varname]]
  if (is.null(v)) stop("Variable not found in netCDF file: ", varname)

  dim_names <- vapply(v$dim, function(d) d$name, character(1))
  i_lon <- match("lon", dim_names)
  i_lat <- match("lat", dim_names)
  i_time <- match("time", dim_names)
  if (any(is.na(c(i_lon, i_lat, i_time)))) {
    stop("Expected dims lon/lat/time for variable ", varname, "; found: ",
         paste(dim_names, collapse = ", "))
  }

  start <- rep(1L, length(dim_names))
  count <- rep(1L, length(dim_names))
  start[i_lon] <- min(lon_idx)
  count[i_lon] <- length(lon_idx)
  start[i_lat] <- min(lat_idx)
  count[i_lat] <- length(lat_idx)
  start[i_time] <- t_index
  count[i_time] <- 1L

  # Let ncdf4 handle scale_factor/add_offset when apply_scale=TRUE.
  x <- ncvar_get(
    nc, varname,
    start = start,
    count = count,
    collapse_degen = FALSE,
    raw_datavals = !apply_scale
  )

  # Reorder to [lat, lon, ...] and drop singleton trailing dims.
  x <- aperm(x, perm = c(i_lat, i_lon, setdiff(seq_along(dim_names), c(i_lat, i_lon))))
  x <- drop(x)
  if (!is.matrix(x)) x <- matrix(x, nrow = length(lat_idx), ncol = length(lon_idx))

  # Apply missing value
  if (!is.na(missing_value)) x[x == missing_value] <- NA

  # Ensure row/col ordering matches lat_keep/lon_keep exactly.
  lat_block <- lat[min(lat_idx):max(lat_idx)]
  lon_block <- lon[min(lon_idx):max(lon_idx)]
  lat_pos <- match(lat_keep, lat_block)
  lon_pos <- match(lon_keep, lon_block)
  x <- x[lat_pos, lon_pos, drop = FALSE]

  return(x)
}

# Vectorize an R x C matrix into length-S vector using manuscript indexing
# (column-major: s = r + (c-1)R).
vec_col_major <- function(mat_rc) {
  as.vector(mat_rc)  # R stores matrices column-major; this matches s=r+(c-1)R
}

# Flag points that lie on land (maps::map.where uses lon in [-180,180]).
is_land_points <- function(lat, lon) {
  if (!requireNamespace("maps", quietly = TRUE)) {
    stop("Package 'maps' is required to exclude land from observation sites.")
  }
  lon180 <- ifelse(lon > 180, lon - 360, lon)
  !is.na(maps::map.where("world", x = lon180, y = lat))
}

# ---------------------------
# 3) Read and subset observations (OISST v2)
# ---------------------------

nc_obs <- nc_open(file_obs)

ym_obs <- get_ym_oisst(nc_obs)
if (!(target_ym %in% ym_obs)) stop("target_ym not found in OISST file. Available range: ",
                                  min(ym_obs), "..", max(ym_obs))

ti_obs <- which(ym_obs == target_ym)[1]

# OISST missing value flag
mv_obs <- ncatt_get(nc_obs, "sst", "missing_value")$value

sst_obs_mat <- extract_window(
  nc = nc_obs,
  varname = "sst",
  t_index = ti_obs,
  lat_keep = lat_keep,
  lon_keep = lon_keep,
  missing_value = mv_obs,
  apply_scale = TRUE
)
nc_close(nc_obs)

# ---------------------------
# 4) Build UQ26-compatible observation objects
# ---------------------------

# Full-lattice vectors
Xobs_full <- vec_col_major(sst_obs_mat)
Yo_full <- Xobs_full

# Build grid metadata consistent with manuscript indexing
# r: 1..R corresponds to lat_keep order (north->south here)
# c: 1..C corresponds to lon_keep order (west->east)

grid_df <- expand.grid(r = 1:R, c = 1:C)
grid_df$s <- grid_df$r + (grid_df$c - 1) * R

grid_df$lat <- rep(lat_keep, times = C)
# For lon, repeat each lon value for all rows
grid_df$lon <- rep(lon_keep, each = R)

# Choose observed sites O from ocean cells only.
land_site <- is_land_points(lat = grid_df$lat, lon = grid_df$lon)
ocean_sites <- grid_df$s[!land_site]
if (N_obs > length(ocean_sites)) {
  stop("N_obs (", N_obs, ") exceeds number of ocean sites (", length(ocean_sites), ").")
}
obs_sites <- sort(sample(ocean_sites, size = N_obs, replace = FALSE))
Yo <- Yo_full[obs_sites]

# Fixed land mask on the current grid.
land_mask_df <- data.frame(
  s = grid_df$s,
  lat = grid_df$lat,
  lon = grid_df$lon,
  is_land = land_site
)

# Build sampling matrix H (sparse) if you want exactly the manuscript form
H <- sparseMatrix(i = 1:N_obs, j = obs_sites, x = 1, dims = c(N_obs, S))

# Bundle for saving
pacific_data <- list(
  meta = list(
    target_ym = target_ym,
    domain = list(lat = lat_keep, lon = lon_keep),
    lattice = list(R = R, C = C, S = S, indexing = "column-major: s = r + (c-1)R"),
    obs_subset = list(N = N_obs, seed = 1)
  ),
  grid = grid_df,
  land_mask = land_mask_df,
  Yo_full = Yo_full,
  Yo = Yo,
  obs_sites = obs_sites,
  H = H
)

# Save outputs
saveRDS(pacific_data, file = file.path(out_dir, "processed", "data.rds"))
saveRDS(pacific_data, file = file.path(out_dir, "processed", "pacific_sst_data.rds"))
write.csv(grid_df, file = file.path(out_dir, "processed", "pacific_sst_grid.csv"), row.names = FALSE)
write.csv(land_mask_df, file = file.path(out_dir, "processed", "pacific_land_mask.csv"), row.names = FALSE)
write.csv(data.frame(s = obs_sites), file = file.path(out_dir, "processed", "pacific_obs_sites.csv"), row.names = FALSE)
write.csv(data.frame(s = 1:S, Yo_full = Yo_full), file = file.path(out_dir, "processed", "pacific_Yo_full.csv"), row.names = FALSE)
write.csv(data.frame(s = obs_sites, Yo = Yo), file = file.path(out_dir, "processed", "pacific_Yo.csv"), row.names = FALSE)

message("Done. Outputs written to: ", normalizePath(out_dir))
