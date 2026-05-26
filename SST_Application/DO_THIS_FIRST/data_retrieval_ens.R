#!/usr/bin/env Rscript

# UQ26: Sea Surface Temperature (SST) ensemble extraction pipeline
#
# This script reads all local E3SM ensemble files:
#   Data/raw/tos_Omon_*.nc
# and extracts the same manuscript-style computer output used by
# DO_THIS_FIRST/data_retrieval_sst.R:
#   - same target month (target_ym)
#   - same lat/lon lattice window
#   - same column-major site indexing s = r + (c-1)R
#
# Outputs:
#   Data/processed/pacific_Yc.csv (ensemble mean, for single-Yc consumers)
#   Data/processed/pacific_Yc_ens_wide.csv
#   Data/processed/pacific_Yc_ens_long.csv
#   Data/processed/pacific_Yc_ens_members.csv
#   Data/processed/pacific_Yc_ensemble_members/*.csv (one pacific_Yc.csv-style file per member)

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
})

out_dir <- "Data"
raw_dir <- file.path(out_dir, "raw")
processed_dir <- file.path(out_dir, "processed")
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

# Defaults aligned with DO_THIS_FIRST/data_retrieval_sst.R.
target_ym <- "2012-01"
lat_keep <- seq(44.5, -44.5, by = -1)
lon_keep <- seq(130.5, 288.5, by = 1)

# Reuse existing processed metadata when available so structure matches exactly.
bundle_file <- file.path(processed_dir, "data.rds")
if (file.exists(bundle_file)) {
  bundle <- readRDS(bundle_file)
  if (!is.null(bundle$meta$target_ym)) target_ym <- as.character(bundle$meta$target_ym)[1]
  if (!is.null(bundle$meta$domain$lat)) lat_keep <- as.numeric(bundle$meta$domain$lat)
  if (!is.null(bundle$meta$domain$lon)) lon_keep <- as.numeric(bundle$meta$domain$lon)
}

R <- length(lat_keep)
C <- length(lon_keep)
S <- R * C

# Map monthly index to YYYY-MM for the 2010-01..2014-12 E3SM files in this project.
get_ym_e3sm <- function() {
  seq(as.Date("2010-01-15"), as.Date("2014-12-15"), by = "1 month") |> format("%Y-%m")
}

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

  x <- ncvar_get(
    nc, varname,
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

member_label <- function(path) {
  b <- basename(path)
  m <- regexec("_r([0-9]+)i[0-9]+p[0-9]+f[0-9]+_", b)
  hit <- regmatches(b, m)[[1]]
  if (length(hit) >= 2L) {
    sprintf("r%02d", as.integer(hit[2]))
  } else {
    tools::file_path_sans_ext(b)
  }
}

member_number <- function(path) {
  b <- basename(path)
  m <- regexec("_r([0-9]+)i[0-9]+p[0-9]+f[0-9]+_", b)
  hit <- regmatches(b, m)[[1]]
  if (length(hit) >= 2L) as.integer(hit[2]) else Inf
}

nc_files <- Sys.glob(file.path(raw_dir, "tos_Omon_*.nc"))
if (length(nc_files) == 0L) {
  stop("No files matched: ", file.path(raw_dir, "tos_Omon_*.nc"))
}

ord <- order(vapply(nc_files, member_number, numeric(1)), basename(nc_files))
nc_files <- nc_files[ord]
members <- vapply(nc_files, member_label, character(1))

yms <- get_ym_e3sm()
if (!(target_ym %in% yms)) {
  stop("target_ym=", target_ym, " is outside expected E3SM range 2010-01..2014-12")
}
t_index <- which(yms == target_ym)[1]

Yc_matrix <- matrix(NA_real_, nrow = S, ncol = length(nc_files))

for (i in seq_along(nc_files)) {
  f <- nc_files[i]
  message(sprintf("[%d/%d] Extracting %s", i, length(nc_files), basename(f)))
  nc <- nc_open(f)
  on.exit(nc_close(nc), add = TRUE)

  mv <- ncatt_get(nc, "tos", "missing_value")$value
  if (is.null(mv) || length(mv) == 0L || is.na(mv[1])) mv <- 1e20

  tos_mat <- extract_window(
    nc = nc,
    varname = "tos",
    t_index = t_index,
    lat_keep = lat_keep,
    lon_keep = lon_keep,
    missing_value = mv[1],
    apply_scale = FALSE
  )

  # Column-major vectorization (s = r + (c-1)R), matching pacific_Yc.csv.
  Yc_matrix[, i] <- as.vector(tos_mat)

  nc_close(nc)
  on.exit(NULL, add = FALSE)
}

colnames(Yc_matrix) <- members

Yc_mean <- rowMeans(Yc_matrix, na.rm = TRUE)

utils::write.csv(
  data.frame(s = seq_len(S), Yc = Yc_mean),
  file = file.path(processed_dir, "pacific_Yc.csv"),
  row.names = FALSE
)

# Full-grid ensemble outputs.
yc_wide_df <- data.frame(s = seq_len(S), Yc_matrix, check.names = FALSE)
utils::write.csv(
  yc_wide_df,
  file = file.path(processed_dir, "pacific_Yc_ens_wide.csv"),
  row.names = FALSE
)

yc_long_df <- data.frame(
  member = rep(members, each = S),
  s = rep(seq_len(S), times = length(members)),
  Yc = as.vector(Yc_matrix),
  stringsAsFactors = FALSE
)
utils::write.csv(
  yc_long_df,
  file = file.path(processed_dir, "pacific_Yc_ens_long.csv"),
  row.names = FALSE
)

utils::write.csv(
  data.frame(member = members, file = basename(nc_files), stringsAsFactors = FALSE),
  file = file.path(processed_dir, "pacific_Yc_ens_members.csv"),
  row.names = FALSE
)

# One pacific_Yc.csv-style file per member.
yc_member_dir <- file.path(processed_dir, "pacific_Yc_ensemble_members")
dir.create(yc_member_dir, showWarnings = FALSE, recursive = TRUE)
for (i in seq_along(members)) {
  utils::write.csv(
    data.frame(s = seq_len(S), Yc = Yc_matrix[, i]),
    file = file.path(yc_member_dir, sprintf("pacific_Yc_%s.csv", members[i])),
    row.names = FALSE
  )
}

message(
  "Done. Processed ensemble outputs written in: ",
  normalizePath(processed_dir, mustWork = FALSE)
)
