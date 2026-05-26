#!/usr/bin/env Rscript

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
  library(Matrix)
})

# Build the ocean-only graph Laplacian from processed SST grid outputs.
# Usage:
#   Rscript DO_THIS_FIRST/compute_ocean_laplacian.R [out_dir] [--no-full-eigen]
#   Rscript DO_THIS_FIRST/compute_ocean_laplacian.R [out_dir] [--force-full-eigen]
#   Rscript DO_THIS_FIRST/compute_ocean_laplacian.R [out_dir] [--max-full-eigen-n=2500]
# Default out_dir:
#   Data

args <- commandArgs(trailingOnly = TRUE)
is_flag <- startsWith(args, "--")
flag_args <- args[is_flag]
positional_args <- args[!is_flag]

out_dir <- if (length(positional_args) >= 1L) positional_args[[1]] else "Data"

compute_full_eigen <- !("--no-full-eigen" %in% flag_args)
force_full_eigen <- "--force-full-eigen" %in% flag_args
max_full_eigen_n <- 2500L

max_n_flags <- grep("^--max-full-eigen-n=", flag_args, value = TRUE)
if (length(max_n_flags) > 1L) {
  stop("At most one --max-full-eigen-n=... flag is allowed.")
}
if (length(max_n_flags) == 1L) {
  max_val <- sub("^--max-full-eigen-n=", "", max_n_flags)
  max_full_eigen_n <- suppressWarnings(as.integer(max_val))
  if (is.na(max_full_eigen_n) || max_full_eigen_n <= 0L) {
    stop("Invalid --max-full-eigen-n value: ", max_val)
  }
}

processed_dir <- file.path(out_dir, "processed")
grid_file <- file.path(processed_dir, "pacific_sst_grid.csv")
land_file <- file.path(processed_dir, "pacific_land_mask.csv")
out_rds <- file.path(processed_dir, "pacific_ocean_graph_laplacian.rds")
out_sites_csv <- file.path(processed_dir, "pacific_ocean_sites_with_local_index.csv")
out_edges_csv <- file.path(processed_dir, "pacific_ocean_edges.csv")
out_eigen_rds <- file.path(processed_dir, "pacific_ocean_laplacian_eigendecomp.rds")

if (!file.exists(grid_file)) stop("Missing grid file: ", grid_file)
if (!file.exists(land_file)) stop("Missing land mask file: ", land_file)

grid_df <- utils::read.csv(grid_file, stringsAsFactors = FALSE)
land_df <- utils::read.csv(land_file, stringsAsFactors = FALSE)

required_grid_cols <- c("s", "r", "c")
required_land_cols <- c("s", "is_land")

if (!all(required_grid_cols %in% names(grid_df))) {
  stop("Grid file must contain columns: ", paste(required_grid_cols, collapse = ", "))
}
if (!all(required_land_cols %in% names(land_df))) {
  stop("Land mask file must contain columns: ", paste(required_land_cols, collapse = ", "))
}

if (anyDuplicated(grid_df$s) > 0L) stop("Duplicate site ids in grid file.")
if (anyDuplicated(land_df$s) > 0L) stop("Duplicate site ids in land mask file.")

site_df <- merge(
  grid_df[, required_grid_cols],
  land_df[, required_land_cols],
  by = "s",
  all.x = TRUE,
  sort = FALSE
)
site_df <- site_df[order(site_df$s), ]

if (anyNA(site_df$is_land)) {
  stop("Land mask does not define all sites in the grid file.")
}

R <- max(site_df$r)
C <- max(site_df$c)
S <- nrow(site_df)

if (R * C != S) {
  stop("Grid shape mismatch: R*C=", R * C, " but number of sites is ", S, ".")
}

expected_s <- site_df$r + (site_df$c - 1L) * R
if (!identical(as.integer(expected_s), as.integer(site_df$s))) {
  stop("Site indexing mismatch. Expected s = r + (c - 1) * R (column-major).")
}

is_land <- as.logical(site_df$is_land)
if (anyNA(is_land)) stop("Could not coerce is_land values to logical.")

ocean_full <- !is_land
ocean_sites <- site_df$s[ocean_full]
n_ocean <- length(ocean_sites)

if (n_ocean == 0L) {
  stop("No ocean cells found in mask.")
}

# Matrix layout follows manuscript indexing: s = r + (c-1)R (column-major).
ocean_rc <- matrix(ocean_full, nrow = R, ncol = C)
s_rc <- matrix(site_df$s, nrow = R, ncol = C)

edge_from <- integer(0)
edge_to <- integer(0)

# Vertical neighbors (N/S): (r, c) <-> (r+1, c), no wrap-around.
if (R > 1L) {
  south_ok <- ocean_rc[1:(R - 1L), , drop = FALSE] & ocean_rc[2:R, , drop = FALSE]
  if (any(south_ok)) {
    edge_from <- c(edge_from, s_rc[1:(R - 1L), , drop = FALSE][south_ok])
    edge_to <- c(edge_to, s_rc[2:R, , drop = FALSE][south_ok])
  }
}

# Horizontal neighbors (E/W): (r, c) <-> (r, c+1), no wrap-around.
if (C > 1L) {
  east_ok <- ocean_rc[, 1:(C - 1L), drop = FALSE] & ocean_rc[, 2:C, drop = FALSE]
  if (any(east_ok)) {
    edge_from <- c(edge_from, s_rc[, 1:(C - 1L), drop = FALSE][east_ok])
    edge_to <- c(edge_to, s_rc[, 2:C, drop = FALSE][east_ok])
  }
}

n_edges <- length(edge_from) # undirected edges counted once

global_to_local <- integer(S)
global_to_local[ocean_sites] <- seq_len(n_ocean)

if (n_edges > 0L) {
  i_local <- global_to_local[edge_from]
  j_local <- global_to_local[edge_to]
  A <- sparseMatrix(
    i = c(i_local, j_local),
    j = c(j_local, i_local),
    x = 1,
    dims = c(n_ocean, n_ocean)
  )
} else {
  A <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n_ocean, n_ocean))
}

degree <- Matrix::rowSums(A)
L <- Diagonal(x = degree) - A

eigen_meta <- list(
  requested = compute_full_eigen,
  computed = FALSE,
  file = NA_character_,
  max_full_eigen_n = as.integer(max_full_eigen_n),
  forced = force_full_eigen
)

if (compute_full_eigen && n_ocean > max_full_eigen_n && !force_full_eigen) {
  warning(
    "Skipping full eigendecomposition for n_ocean=", n_ocean,
    " because n_ocean > max_full_eigen_n (", max_full_eigen_n, "). ",
    "Use --force-full-eigen (or raise --max-full-eigen-n) to attempt it anyway."
  )
  compute_full_eigen <- FALSE
}

if (compute_full_eigen) {
  message("Computing full eigendecomposition: L = U %*% Lambda %*% t(U)")
  L_dense <- as.matrix(L)
  eig <- eigen(L_dense, symmetric = TRUE)

  # Reconstruct with elementwise scaling of t(U) columns by eigenvalues.
  recon <- eig$vectors %*% (eig$values * t(eig$vectors))
  max_abs_reconstruction_error <- max(abs(L_dense - recon))

  eigen_obj <- list(
    U = eig$vectors,
    Lambda = Diagonal(x = eig$values),
    eigenvalues = eig$values,
    max_abs_reconstruction_error = max_abs_reconstruction_error
  )
  saveRDS(eigen_obj, out_eigen_rds)

  eigen_meta$computed <- TRUE
  eigen_meta$file <- out_eigen_rds
  eigen_meta$max_abs_reconstruction_error <- max_abs_reconstruction_error

  rm(L_dense, eig, recon, eigen_obj)
  invisible(gc())
}

ocean_index_df <- data.frame(
  local_index = seq_len(n_ocean),
  s = ocean_sites,
  r = site_df$r[ocean_full],
  c = site_df$c[ocean_full]
)

edge_df <- data.frame(s_from = edge_from, s_to = edge_to)

res <- list(
  meta = list(
    R = R,
    C = C,
    S = S,
    n_ocean = n_ocean,
    n_land = S - n_ocean,
    n_edges_undirected = n_edges,
    indexing = "column-major: s = r + (c - 1) * R",
    neighborhood = "N/S/E/W, no wrap-around"
  ),
  ocean_sites = ocean_sites,
  ocean_index = ocean_index_df,
  adjacency = A,
  degree = as.numeric(degree),
  laplacian = L,
  eigendecomposition = eigen_meta
)

saveRDS(res, out_rds)
utils::write.csv(ocean_index_df, out_sites_csv, row.names = FALSE)
utils::write.csv(edge_df, out_edges_csv, row.names = FALSE)

message(
  "Wrote ocean-graph Laplacian to: ", normalizePath(out_rds), "\n",
  "Ocean cells: ", n_ocean, " / ", S, "\n",
  "Undirected edges: ", n_edges, "\n",
  "Full eigendecomposition computed: ", eigen_meta$computed,
  if (isTRUE(eigen_meta$computed)) {
    paste0("\nEigendecomposition file: ", normalizePath(out_eigen_rds))
  } else {
    ""
  }
)
