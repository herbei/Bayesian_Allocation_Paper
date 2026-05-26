# Build SST problem objects for the spectral MCMC workflow.
#
# Returns:
#   list(
#     static = ...   # object consumed by modules/mcmc_*.R
#     setup_snapshot = ...
#     meta = ...
#   )
#
# The returned "static" object uses the same field names expected by the
# original spectral2 module functions (S, U, lambda, Yc, Yo, calH, ...).
build_sst_problem_data <- function(data_dir = "Data", K_modes_target = 500L) {
  processed_dir <- file.path(data_dir, "processed")

  grid_file <- file.path(processed_dir, "pacific_sst_grid.csv")
  land_file <- file.path(processed_dir, "pacific_land_mask.csv")
  yc_file <- file.path(processed_dir, "pacific_Yc.csv")
  yo_file <- file.path(processed_dir, "pacific_Yo.csv")
  lap_file <- file.path(processed_dir, "pacific_ocean_graph_laplacian.rds")
  eig_file <- file.path(processed_dir, "pacific_ocean_laplacian_eigendecomp.rds")

  required_files <- c(grid_file, land_file, yc_file, yo_file, lap_file, eig_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop("Missing required SST inputs:\n", paste0("  - ", missing_files, collapse = "\n"))
  }

  grid_df <- utils::read.csv(grid_file, stringsAsFactors = FALSE)
  land_df <- utils::read.csv(land_file, stringsAsFactors = FALSE)
  yc_df <- utils::read.csv(yc_file, stringsAsFactors = FALSE)
  yo_df <- utils::read.csv(yo_file, stringsAsFactors = FALSE)
  lap_obj <- readRDS(lap_file)
  eig_obj <- readRDS(eig_file)

  required_grid_cols <- c("s", "r", "c", "lat", "lon")
  required_land_cols <- c("s", "is_land")
  required_yc_cols <- c("s", "Yc")
  required_yo_cols <- c("s", "Yo")

  if (!all(required_grid_cols %in% names(grid_df))) {
    stop("Grid file must contain: ", paste(required_grid_cols, collapse = ", "))
  }
  if (!all(required_land_cols %in% names(land_df))) {
    stop("Land-mask file must contain: ", paste(required_land_cols, collapse = ", "))
  }
  if (!all(required_yc_cols %in% names(yc_df))) {
    stop("Yc file must contain: ", paste(required_yc_cols, collapse = ", "))
  }
  if (!all(required_yo_cols %in% names(yo_df))) {
    stop("Yo file must contain: ", paste(required_yo_cols, collapse = ", "))
  }

  grid_df <- grid_df[order(grid_df$s), required_grid_cols]
  land_df <- land_df[order(land_df$s), required_land_cols]
  yc_df <- yc_df[order(yc_df$s), required_yc_cols]
  yo_df <- yo_df[order(yo_df$s), required_yo_cols]

  S_full <- nrow(grid_df)
  R_full <- length(unique(grid_df$r))
  C_full <- length(unique(grid_df$c))

  if (R_full * C_full != S_full) {
    stop("Grid dimensions inconsistent: R*C != number of grid cells.")
  }
  if (!identical(as.integer(grid_df$s), seq_len(S_full))) {
    stop("Expected grid ids s = 1..S.")
  }

  expected_s <- grid_df$r + (grid_df$c - 1L) * R_full
  if (!identical(as.integer(expected_s), as.integer(grid_df$s))) {
    stop("Grid indexing mismatch: expected s = r + (c - 1) * R.")
  }
  if (!identical(as.integer(land_df$s), as.integer(grid_df$s))) {
    stop("Land-mask ids do not match grid ids.")
  }
  if (!identical(as.integer(yc_df$s), as.integer(grid_df$s))) {
    stop("Yc ids do not match grid ids.")
  }

  ocean_sites <- as.integer(lap_obj$ocean_sites)
  n_ocean <- length(ocean_sites)
  if (n_ocean == 0L) stop("No ocean sites found in ocean graph object.")

  is_land <- as.logical(land_df$is_land)
  if (anyNA(is_land)) stop("Could not coerce land mask values to logical.")
  if (sum(!is_land) != n_ocean) {
    stop("Land-mask ocean count does not match ocean graph ocean count.")
  }

  global_to_local <- integer(S_full)
  global_to_local[ocean_sites] <- seq_len(n_ocean)

  # Full model output restricted to ocean cells.
  Yc_ocean <- as.numeric(yc_df$Yc)[ocean_sites]
  bad_yc <- which(!is.finite(Yc_ocean))
  yc_fill <- NA_real_
  if (length(bad_yc) > 0L) {
    yc_fill <- stats::median(Yc_ocean[is.finite(Yc_ocean)])
    Yc_ocean[bad_yc] <- yc_fill
    warning(sprintf(
      "Imputed %d non-finite ocean Yc values with ocean median %.6f.",
      length(bad_yc), yc_fill
    ))
  }

  # Observations are sampled on ocean sites in DO_THIS_FIRST/data_retrieval_sst.R.
  obs_sites_global <- as.integer(yo_df$s)
  obs_sites_local <- global_to_local[obs_sites_global]
  if (any(obs_sites_local <= 0L)) {
    stop("Found observed sites that are not ocean sites.")
  }

  Yo_obs <- as.numeric(yo_df$Yo)
  if (any(!is.finite(Yo_obs))) {
    stop("Yo contains non-finite values; please clean observations first.")
  }
  N_obs <- length(Yo_obs)

  # Use a dense observation matrix so runtime does not require the Matrix package.
  # (N_obs x n_ocean = 1000 x 12554 is manageable for current settings.)
  calH <- matrix(0, nrow = N_obs, ncol = n_ocean)
  calH[cbind(seq_len(N_obs), obs_sites_local)] <- 1

  if (is.null(eig_obj$U) || is.null(eig_obj$eigenvalues)) {
    stop("Eigendecomposition object must contain U and eigenvalues.")
  }

  lambda_all <- as.numeric(eig_obj$eigenvalues)
  if (length(lambda_all) != n_ocean) {
    stop("Eigenvalue count does not match number of ocean sites.")
  }
  if (!all(dim(eig_obj$U) == c(n_ocean, n_ocean))) {
    stop("Eigenvector matrix dimension mismatch.")
  }

  ord <- order(lambda_all, decreasing = FALSE)
  K_modes <- max(2L, min(as.integer(K_modes_target), n_ocean))
  mode_indices <- ord[seq_len(K_modes)]

  lambda_modes <- lambda_all[mode_indices]
  U_modes <- eig_obj$U[, mode_indices, drop = FALSE]

  rm(eig_obj)
  invisible(gc())

  # Match the field names expected by modules/mcmc_*.R.
  static <- list(
    S = n_ocean,
    R = R_full,
    C = C_full,
    U = U_modes,
    lambda = lambda_modes,
    Yc = as.numeric(Yc_ocean),
    Yo = as.numeric(Yo_obs),
    Yo_obs = as.numeric(Yo_obs),
    calH = calH,
    N = N_obs,
    U_obs = U_modes[obs_sites_local, , drop = FALSE],
    Zc_all = as.numeric(crossprod(U_modes, Yc_ocean))
  )

  setup_snapshot <- list(
    data_dir = normalizePath(data_dir, mustWork = FALSE),
    grid_dims = list(R = R_full, C = C_full, S_full = S_full),
    n_ocean = n_ocean,
    obs_count = N_obs,
    ocean_sites_global = ocean_sites,
    obs_sites_global = obs_sites_global,
    obs_sites_local = obs_sites_local,
    mode_indices_in_full_basis = mode_indices,
    imputed_nonfinite_Yc_count = length(bad_yc),
    imputed_Yc_value = yc_fill
  )

  meta <- list(
    K_modes = K_modes,
    n_ocean = n_ocean,
    N_obs = N_obs,
    R = R_full,
    C = C_full,
    S_full = S_full
  )

  list(static = static, setup_snapshot = setup_snapshot, meta = meta)
}
