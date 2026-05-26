rm(list = ls())

##########
# 1) Environment and module loading
##########

source("modules/sst_data_setup.R")
source("modules/mcmc_utils.R")
source("modules/mcmc_state.R")
source("modules/mcmc_updates.R")
source("modules/mcmc_storage.R")
source("modules/mcmc_run.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

member_number <- function(path) {
  b <- basename(path)
  patterns <- c(
    "_r([0-9]+)i[0-9]+p[0-9]+f[0-9]+_",
    "pacific_Yc_r([0-9]+)",
    "_r([0-9]+)\\."
  )
  for (pat in patterns) {
    m <- regexec(pat, b, perl = TRUE)
    hit <- regmatches(b, m)[[1]]
    if (length(hit) >= 2L) return(as.integer(hit[2]))
  }
  Inf
}

member_label <- function(path) {
  idx <- member_number(path)
  if (is.finite(idx)) return(sprintf("r%02d", as.integer(idx)))
  tools::file_path_sans_ext(basename(path))
}

# Build SST problem objects for the spectral MCMC workflow (ensemble version).
#
# Returns:
#   list(
#     static = ...   # object consumed by run_mcmc_ens and update functions
#     setup_snapshot = ...
#     meta = ...
#   )
#
# Ensemble inputs are read from:
#   Data/processed/pacific_Yc_ensemble_members/pacific_Yc_*.csv
build_sst_problem_data_ens <- function(data_dir = "Data", K_modes_target = 500L) {
  processed_dir <- file.path(data_dir, "processed")

  grid_file <- file.path(processed_dir, "pacific_sst_grid.csv")
  land_file <- file.path(processed_dir, "pacific_land_mask.csv")
  yo_file <- file.path(processed_dir, "pacific_Yo.csv")
  lap_file <- file.path(processed_dir, "pacific_ocean_graph_laplacian.rds")
  eig_file <- file.path(processed_dir, "pacific_ocean_laplacian_eigendecomp.rds")
  yc_dir <- file.path(processed_dir, "pacific_Yc_ensemble_members")
  yc_files <- Sys.glob(file.path(yc_dir, "pacific_Yc_*.csv"))

  required_files <- c(grid_file, land_file, yo_file, lap_file, eig_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop("Missing required SST inputs:\n", paste0("  - ", missing_files, collapse = "\n"))
  }
  if (length(yc_files) == 0L) {
    stop("No ensemble member files found in ", yc_dir, " (expected pacific_Yc_*.csv).")
  }

  ord_files <- order(vapply(yc_files, member_number, numeric(1)), basename(yc_files))
  yc_files <- yc_files[ord_files]
  yc_members <- vapply(yc_files, member_label, character(1))
  M_ens <- length(yc_files)

  grid_df <- utils::read.csv(grid_file, stringsAsFactors = FALSE)
  land_df <- utils::read.csv(land_file, stringsAsFactors = FALSE)
  yo_df <- utils::read.csv(yo_file, stringsAsFactors = FALSE)
  lap_obj <- readRDS(lap_file)
  eig_obj <- readRDS(eig_file)

  required_grid_cols <- c("s", "r", "c", "lat", "lon")
  required_land_cols <- c("s", "is_land")
  required_yo_cols <- c("s", "Yo")

  if (!all(required_grid_cols %in% names(grid_df))) {
    stop("Grid file must contain: ", paste(required_grid_cols, collapse = ", "))
  }
  if (!all(required_land_cols %in% names(land_df))) {
    stop("Land-mask file must contain: ", paste(required_land_cols, collapse = ", "))
  }
  if (!all(required_yo_cols %in% names(yo_df))) {
    stop("Yo file must contain: ", paste(required_yo_cols, collapse = ", "))
  }

  grid_df <- grid_df[order(grid_df$s), required_grid_cols]
  land_df <- land_df[order(land_df$s), required_land_cols]
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

  # Read full-grid Yc member fields and align by site id.
  Yc_full_matrix <- matrix(NA_real_, nrow = S_full, ncol = M_ens)
  colnames(Yc_full_matrix) <- yc_members

  for (j in seq_along(yc_files)) {
    yc_df <- utils::read.csv(yc_files[j], stringsAsFactors = FALSE)
    if (!all(c("s", "Yc") %in% names(yc_df))) {
      stop("Yc ensemble file must contain columns s,Yc: ", yc_files[j])
    }
    yc_df <- yc_df[order(yc_df$s), c("s", "Yc")]
    if (!identical(as.integer(yc_df$s), as.integer(grid_df$s))) {
      stop("Yc ensemble ids do not match grid ids for file: ", yc_files[j])
    }
    Yc_full_matrix[, j] <- as.numeric(yc_df$Yc)
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

  # Restrict ensemble model output to ocean cells and impute non-finite values per member.
  Yc_ocean_matrix <- Yc_full_matrix[ocean_sites, , drop = FALSE]
  rm(Yc_full_matrix)
  invisible(gc())

  imputed_nonfinite_per_member <- integer(M_ens)
  imputed_fill_per_member <- rep(NA_real_, M_ens)
  for (j in seq_len(M_ens)) {
    yj <- as.numeric(Yc_ocean_matrix[, j])
    bad_j <- which(!is.finite(yj))
    imputed_nonfinite_per_member[j] <- length(bad_j)
    if (length(bad_j) > 0L) {
      fill_j <- stats::median(yj[is.finite(yj)])
      if (!is.finite(fill_j)) {
        fill_j <- 0
        warning(sprintf(
          "Member %s had no finite ocean Yc values; imputed all with 0.",
          yc_members[j]
        ))
      }
      yj[bad_j] <- fill_j
      imputed_fill_per_member[j] <- fill_j
      Yc_ocean_matrix[, j] <- yj
    }
  }

  if (any(imputed_nonfinite_per_member > 0L)) {
    nz <- which(imputed_nonfinite_per_member > 0L)
    warning(sprintf(
      "Imputed %d non-finite ocean Yc values across %d ensemble members.",
      sum(imputed_nonfinite_per_member[nz]),
      length(nz)
    ))
  }

  Yc_bar_ocean <- rowMeans(Yc_ocean_matrix)
  Wc <- sum((Yc_ocean_matrix - Yc_bar_ocean)^2)

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

  # Dense observation matrix, same style as single-output setup.
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

  # Match the field names expected by existing mcmc_* style code.
  U_obs <- U_modes[obs_sites_local, , drop = FALSE]
  G_obs <- crossprod(U_obs)
  Yo_proj_all <- as.numeric(crossprod(U_obs, Yo_obs))

  static <- list(
    S = n_ocean,
    R = R_full,
    C = C_full,
    U = U_modes,
    lambda = lambda_modes,
    Yc = as.numeric(Yc_bar_ocean),    # Ensemble mean ybar^c
    Yc_bar = as.numeric(Yc_bar_ocean),
    M_ens = as.integer(M_ens),
    Wc = as.numeric(Wc),
    Yo = as.numeric(Yo_obs),
    Yo_obs = as.numeric(Yo_obs),
    calH = calH,
    N = N_obs,
    U_obs = U_obs,
    G_obs = G_obs,
    U_obs_col_sq = as.numeric(diag(G_obs)),
    Yo_proj_all = Yo_proj_all,
    Zc_all = as.numeric(crossprod(U_modes, Yc_bar_ocean))
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
    ensemble_count = M_ens,
    ensemble_member_labels = yc_members,
    ensemble_member_files = basename(yc_files),
    Wc = as.numeric(Wc),
    imputed_nonfinite_per_member = imputed_nonfinite_per_member,
    imputed_fill_per_member = imputed_fill_per_member
  )

  meta <- list(
    K_modes = K_modes,
    n_ocean = n_ocean,
    N_obs = N_obs,
    R = R_full,
    C = C_full,
    S_full = S_full,
    M_ens = M_ens
  )

  list(static = static, setup_snapshot = setup_snapshot, meta = meta)
}

# Ensemble generalization of update_B_local (Section 6.4): replace sigma_c^-2 by M*sigma_c^-2.
update_B_local_ens <- function(B_vec, sig2c_curr, sig2o_curr, static, params) {
  inv_sig2c_eff <- static$M_ens / sig2c_curr
  inv_sig2o <- 1 / sig2o_curr
  inv_sig2c_eff2 <- inv_sig2c_eff^2
  eps_num <- 1e-12
  G_obs <- static$G_obs
  col_sq <- static$U_obs_col_sq
  Yo_proj_all <- static$Yo_proj_all

  I_X_work <- as.integer(which(B_vec == 1L))
  Z_X_work <- static$Zc_all[I_X_work]
  h_X_work <- inv_sig2c_eff * Z_X_work + inv_sig2o * Yo_proj_all[I_X_work]

  if (length(I_X_work) > 0L) {
    P_X_work <- G_obs[I_X_work, I_X_work, drop = FALSE] * inv_sig2o
    D_X_work <- params$delta2 * static$lambda[I_X_work] + params$delta02
    diag(P_X_work) <- diag(P_X_work) + D_X_work + inv_sig2c_eff
    # One factorization at the start of the sweep; then update inverse locally.
    P_X_inv_work <- chol2inv(chol(P_X_work))
    P_X_inv_work <- 0.5 * (P_X_inv_work + t(P_X_inv_work))
  } else {
    P_X_work <- matrix(0, nrow = 0, ncol = 0)
    P_X_inv_work <- matrix(0, nrow = 0, ncol = 0)
  }

  for (k in params$calK_0) {
    pk_curr <- as.numeric(params$p_k[as.character(k)])
    pk_curr <- min(max(pk_curr, eps_num), 1 - eps_num)

    z_k <- static$Zc_all[k]
    d_X_k <- params$delta2 * static$lambda[k] + params$delta02
    d_beta_k <- params$delta1 * static$lambda[k] + params$delta01

    beta_term_k <- 0.5 * (
      log(d_beta_k) - log(d_beta_k + inv_sig2c_eff) +
      (inv_sig2c_eff2 * z_k^2) / (d_beta_k + inv_sig2c_eff)
    )

    a_k_sq <- col_sq[k]
    y_a_k <- Yo_proj_all[k]

    if (B_vec[k] == 1L) {
      idx_k <- match(k, I_X_work)
      d_X_old <- length(I_X_work)

      if (d_X_old > 1L) {
        alpha_k <- as.numeric(P_X_inv_work[idx_k, idx_k])
        if (!is.finite(alpha_k) || alpha_k <= eps_num) {
          # Rare numerical fallback: refresh inverse from current matrix.
          P_X_inv_work <- chol2inv(chol(P_X_work))
          P_X_inv_work <- 0.5 * (P_X_inv_work + t(P_X_inv_work))
          alpha_k <- as.numeric(P_X_inv_work[idx_k, idx_k])
        }

        s_k <- 1 / alpha_k
        xk <- sum(P_X_inv_work[idx_k, ] * h_X_work)
        g_k <- xk / alpha_k
      } else {
        s_k <- as.numeric(P_X_work[1, 1])
        g_k <- as.numeric(h_X_work[1])
      }

      s_k <- max(s_k, eps_num)
      delta_X_10 <- 0.5 * (-log(d_X_k) + log(s_k) - (g_k^2) / s_k)
      delta_Psi_10 <- log((1 - pk_curr) / pk_curr) + beta_term_k + delta_X_10

      prob_one <- inv_logit(-delta_Psi_10)
      B_new <- as.integer(rbinom(1, size = 1, prob = prob_one))

      if (B_new == 0L) {
        B_vec[k] <- 0L

        if (d_X_old > 1L) {
          keep_idx <- seq_len(d_X_old)[-idx_k]
          alpha_k <- as.numeric(P_X_inv_work[idx_k, idx_k])
          p_col <- P_X_inv_work[keep_idx, idx_k, drop = FALSE]

          # Downdate inverse using block-inverse identity.
          P_X_inv_work <- P_X_inv_work[keep_idx, keep_idx, drop = FALSE] -
            (p_col %*% t(p_col)) / alpha_k
          P_X_inv_work <- 0.5 * (P_X_inv_work + t(P_X_inv_work))

          P_X_work <- P_X_work[keep_idx, keep_idx, drop = FALSE]
          I_X_work <- I_X_work[keep_idx]
          Z_X_work <- Z_X_work[keep_idx]
          h_X_work <- h_X_work[keep_idx]
        } else {
          P_X_work <- matrix(0, nrow = 0, ncol = 0)
          P_X_inv_work <- matrix(0, nrow = 0, ncol = 0)
          I_X_work <- integer(0)
          Z_X_work <- numeric(0)
          h_X_work <- numeric(0)
        }
      }
    } else {
      d_X_old <- length(I_X_work)
      a_add <- d_X_k + inv_sig2c_eff + inv_sig2o * a_k_sq
      h_add <- inv_sig2c_eff * z_k + inv_sig2o * y_a_k

      if (d_X_old > 0) {
        b_add <- inv_sig2o * G_obs[I_X_work, k]
        C_inv_b <- as.numeric(P_X_inv_work %*% b_add)
        C_inv_hr <- as.numeric(P_X_inv_work %*% h_X_work)

        s_k <- as.numeric(a_add - sum(b_add * C_inv_b))
        g_k <- as.numeric(h_add - sum(b_add * C_inv_hr))
      } else {
        b_add <- numeric(0)
        C_inv_b <- numeric(0)
        s_k <- as.numeric(a_add)
        g_k <- as.numeric(h_add)
      }

      s_k <- max(s_k, eps_num)
      delta_X_01 <- 0.5 * (log(d_X_k) - log(s_k) + (g_k^2) / s_k)
      delta_Psi_01 <- log(pk_curr / (1 - pk_curr)) - beta_term_k + delta_X_01

      prob_one <- inv_logit(delta_Psi_01)
      B_new <- as.integer(rbinom(1, size = 1, prob = prob_one))

      if (B_new == 1L) {
        B_vec[k] <- 1L

        I_X_work <- c(I_X_work, k)
        Z_X_work <- c(Z_X_work, z_k)
        h_X_work <- c(h_X_work, h_add)

        if (d_X_old > 0) {
          # Rank-1 block update of inverse for appended mode.
          P11_new <- P_X_inv_work + tcrossprod(C_inv_b) / s_k
          P12_new <- -C_inv_b / s_k
          P_X_inv_work <- rbind(
            cbind(P11_new, P12_new),
            c(P12_new, 1 / s_k)
          )
          P_X_inv_work <- 0.5 * (P_X_inv_work + t(P_X_inv_work))

          P_X_work <- rbind(
            cbind(P_X_work, b_add),
            c(b_add, a_add)
          )
        } else {
          P_X_work <- matrix(a_add, nrow = 1, ncol = 1)
          P_X_inv_work <- matrix(1 / s_k, nrow = 1, ncol = 1)
        }
      }
    }
  }

  B_vec
}

build_state_from_B_ens <- function(B_vec, static, params) {
  n_modes <- ncol(static$U)
  if (length(B_vec) != n_modes) {
    stop(sprintf("Length of B_vec (%d) must match number of modes (%d).", length(B_vec), n_modes))
  }

  I_beta_loc <- as.integer(which(B_vec == 0L))
  I_X_loc <- as.integer(setdiff(seq_len(n_modes), I_beta_loc))

  d_beta_loc <- length(I_beta_loc)
  d_X_loc <- length(I_X_loc)

  U_beta_loc <- static$U[, I_beta_loc, drop = FALSE]
  U_X_loc <- static$U[, I_X_loc, drop = FALSE]

  D_beta_loc <- params$delta1 * static$lambda[I_beta_loc] + params$delta01
  D_X_loc <- params$delta2 * static$lambda[I_X_loc] + params$delta02

  # Reuse precomputed Laplacian coefficients of Ybar^c instead of recomputing crossprod each sweep.
  Zc_beta_loc <- if (d_beta_loc > 0) as.numeric(static$Zc_all[I_beta_loc]) else numeric(0)
  Zc_X_loc <- if (d_X_loc > 0) as.numeric(static$Zc_all[I_X_loc]) else numeric(0)

  list(
    I_beta = I_beta_loc,
    I_X = I_X_loc,
    d_beta = d_beta_loc,
    d_X = d_X_loc,
    U_beta = U_beta_loc,
    U_X = U_X_loc,
    D_beta = D_beta_loc,
    D_X = D_X_loc,
    Zc_beta = Zc_beta_loc,
    Zc_X = Zc_X_loc
  )
}

sample_b_ens <- function(state, sig2c_curr, static) {
  d_beta <- state$d_beta
  inv_sig2c_eff <- static$M_ens / sig2c_curr

  if (d_beta > 0) {
    Pb_diag <- state$D_beta + inv_sig2c_eff
    rhs_b <- inv_sig2c_eff * state$Zc_beta
    m_b <- rhs_b / Pb_diag
    b <- as.numeric(m_b + rnorm(d_beta) / sqrt(Pb_diag))
  } else {
    Pb_diag <- numeric(0)
    m_b <- numeric(0)
    b <- numeric(0)
  }

  list(b = b, m_b = m_b, Pb_diag = Pb_diag)
}

sample_x_ens <- function(state, sig2c_curr, sig2o_curr, static) {
  d_X <- state$d_X
  # calH is a row-selector matrix; U_obs already stores the selected rows.
  A_X <- static$U_obs[, state$I_X, drop = FALSE]
  inv_sig2c_eff <- static$M_ens / sig2c_curr

  if (d_X > 0) {
    Px <- static$G_obs[state$I_X, state$I_X, drop = FALSE] / sig2o_curr
    diag(Px) <- diag(Px) + state$D_X + inv_sig2c_eff

    rhs_x <- inv_sig2c_eff * state$Zc_X + static$Yo_proj_all[state$I_X] / sig2o_curr
    Px_chol <- chol(Px)
    m_x <- as.numeric(backsolve(Px_chol, forwardsolve(t(Px_chol), rhs_x)))
    x <- as.numeric(m_x + backsolve(Px_chol, rnorm(d_X)))
  } else {
    Px_chol <- matrix(0, nrow = 0, ncol = 0)
    m_x <- numeric(0)
    x <- numeric(0)
  }

  list(x = x, m_x = m_x, A_X = A_X, Px_chol = Px_chol)
}

update_sigmas_ens <- function(state, b, x, A_X, static, params) {
  X_mean_obs <- as.numeric(A_X %*% x)
  sse_o <- sum((static$Yo_obs - X_mean_obs)^2)
  sig2o_new <- rinvgamma1(params$gamma_o + static$N / 2, params$xi_o + sse_o / 2)

  beta_vec <- if (state$d_beta > 0) as.numeric(state$U_beta %*% b) else rep(0, static$S)
  X_vec <- as.numeric(state$U_X %*% x)

  sse_c_mean <- sum((static$Yc_bar - beta_vec - X_vec)^2)
  sse_c <- static$Wc + static$M_ens * sse_c_mean
  sig2c_new <- rinvgamma1(
    params$gamma_c + (static$M_ens * static$S) / 2,
    params$xi_c + sse_c / 2
  )

  list(sig2o = sig2o_new, sig2c = sig2c_new, beta_vec = beta_vec, X_vec = X_vec)
}

compute_log_mcmc_summaries_ens <- function(state, x_curr, b_curr, A_X, static, sig2o_curr, sig2c_curr) {
  two_pi <- 2 * pi

  X_mean_obs <- as.numeric(A_X %*% x_curr)
  resid_o <- static$Yo_obs - X_mean_obs
  logp_Yo <- 0.5 * static$N * log(two_pi * sig2o_curr) + 0.5 * sum(resid_o^2) / sig2o_curr

  sig2c_eff <- sig2c_curr / static$M_ens

  if (state$d_X > 0) {
    resid_ZXc <- state$Zc_X - x_curr
    logp_ZXc <- 0.5 * state$d_X * log(two_pi * sig2c_eff) + 0.5 * sum(resid_ZXc^2) / sig2c_eff
  } else {
    logp_ZXc <- 0
  }

  if (state$d_beta > 0) {
    resid_Zbetac <- state$Zc_beta - b_curr
    logp_Zbetac <- 0.5 * state$d_beta * log(two_pi * sig2c_eff) + 0.5 * sum(resid_Zbetac^2) / sig2c_eff
  } else {
    logp_Zbetac <- 0
  }

  list(
    logp_Yo = logp_Yo,
    logp_ZXc = logp_ZXc,
    logp_Zbetac = logp_Zbetac
  )
}

run_mcmc_ens <- function(B_init, sig2o_init, sig2c_init, n_iter, burn_in,
                         progress_every, MCstore, static, params,
                         update_sigma2yo = TRUE,
                         update_sigma2yc = TRUE,
                         update_sigma2_after = 15000L,
                         save_rdata_path = NULL,
                         setup_snapshot = NULL) {
  if (length(update_sigma2_after) != 1L || is.na(update_sigma2_after)) {
    stop("update_sigma2_after must be a single non-missing integer >= 0.")
  }
  update_sigma2_after <- max(0L, as.integer(update_sigma2_after))
  update_sigma2yo <- isTRUE(update_sigma2yo)
  update_sigma2yc <- isTRUE(update_sigma2yc)

  store <- init_storage(n_iter, MCstore, static, params)

  B_curr <- B_init
  sig2o_curr <- sig2o_init
  sig2c_curr <- sig2c_init

  state_curr <- build_state_from_B_ens(B_curr, static, params)
  b_curr <- rep(0, state_curr$d_beta)
  x_curr <- rep(0, state_curr$d_X)

  for (iter in seq_len(n_iter)) {
    # Step 1: local one-mode collapsed update for B (ensemble precision M/sig2c).
    B_curr <- update_B_local_ens(B_curr, sig2c_curr, sig2o_curr, static, params)

    # Step 2: rebuild state under updated B.
    state_curr <- build_state_from_B_ens(B_curr, static, params)

    # Step 3: sample b | B, sig2c, Ybar^c.
    b_res <- sample_b_ens(state_curr, sig2c_curr, static)
    b_curr <- b_res$b

    # Step 4: sample x | B, sig2c, sig2o, Ybar^c, Yo.
    x_res <- sample_x_ens(state_curr, sig2c_curr, sig2o_curr, static)
    x_curr <- x_res$x

    # Step 5/6: sample sig2o and/or sig2c after the requested warmup period.
    do_update_sigma2yo <- update_sigma2yo && (iter > update_sigma2_after)
    do_update_sigma2yc <- update_sigma2yc && (iter > update_sigma2_after)
    if (do_update_sigma2yo || do_update_sigma2yc) {
      sig_res <- update_sigmas_ens(state_curr, b_curr, x_curr, x_res$A_X, static, params)
      if (do_update_sigma2yo) {
        sig2o_curr <- sig_res$sig2o
      }
      if (do_update_sigma2yc) {
        sig2c_curr <- sig_res$sig2c
      }
      beta_vec_curr <- sig_res$beta_vec
      X_vec_curr <- sig_res$X_vec
    } else {
      beta_vec_curr <- if (state_curr$d_beta > 0) as.numeric(state_curr$U_beta %*% b_curr) else rep(0, static$S)
      X_vec_curr <- as.numeric(state_curr$U_X %*% x_curr)
    }

    log_summ <- compute_log_mcmc_summaries_ens(
      state = state_curr,
      x_curr = x_curr,
      b_curr = b_curr,
      A_X = x_res$A_X,
      static = static,
      sig2o_curr = sig2o_curr,
      sig2c_curr = sig2c_curr
    )

    # Step 7: store if requested by MCstore.
    store <- store_draw(
      store = store,
      iter = iter,
      B_vec = B_curr,
      state = state_curr,
      sig2o_curr = sig2o_curr,
      sig2c_curr = sig2c_curr,
      logp_Yo_curr = log_summ$logp_Yo,
      logp_ZXc_curr = log_summ$logp_ZXc,
      logp_Zbetac_curr = log_summ$logp_Zbetac,
      beta_vec = beta_vec_curr,
      X_vec = X_vec_curr,
      b_vec = b_curr,
      x_vec = x_curr,
      params = params
    )

    if (iter %% progress_every == 0 || iter == 1 || iter == n_iter) {
      cat(sprintf(
        "iter %d/%d: d_beta=%d, sig2o=%.4f, sig2c=%.4f, nlogYo=%.2f, nlogZXc=%.2f, nlogZbetac=%.2f\n",
        iter, n_iter, state_curr$d_beta, sig2o_curr, sig2c_curr,
        log_summ$logp_Yo, log_summ$logp_ZXc, log_summ$logp_Zbetac
      ))
    }
  }

  keep_idx <- which(store$iter_chain > burn_in)
  if (length(keep_idx) == 0L) {
    keep_idx <- store$n_store
  }

  beta_post_mean <- rowMeans(store$beta_chain[, keep_idx, drop = FALSE])
  X_post_mean <- rowMeans(store$X_chain[, keep_idx, drop = FALSE])
  phi_post_mean <- beta_post_mean + X_post_mean

  mcmc_output <- list(
    store = store,
    keep_idx = keep_idx,
    beta_post_mean = beta_post_mean,
    X_post_mean = X_post_mean,
    phi_post_mean = phi_post_mean
  )

  if (!is.null(save_rdata_path) && nzchar(save_rdata_path)) {
    save_dir <- dirname(save_rdata_path)
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    }

    mcmc_artifact <- list(
      meta = list(
        saved_at = as.character(Sys.time()),
        source = "run_mcmc_ens",
        format_version = 1L
      ),
      run_config = list(
        n_iter = n_iter,
        burn_in = burn_in,
        progress_every = progress_every,
        MCstore = MCstore,
        update_sigma2yo = update_sigma2yo,
        update_sigma2yc = update_sigma2yc,
        update_sigma2_after = update_sigma2_after,
        sig2o_init = sig2o_init,
        sig2c_init = sig2c_init
      ),
      initial_state = list(
        B_init = B_init
      ),
      static = static,
      params = params,
      setup = setup_snapshot,
      output = mcmc_output
    )

    save(mcmc_artifact, file = save_rdata_path)
    mcmc_output$save_rdata_path <- normalizePath(save_rdata_path, mustWork = FALSE)
    cat(sprintf("Saved MCMC artifact to %s\n", mcmc_output$save_rdata_path))
  }

  mcmc_output
}

##########
# 2) Data loading and graph-basis setup (SST ocean graph, ensemble Yc)
##########

data_dir <- getOption("sst_spectral_ens_data_dir", getOption("sst_spectral_data_dir", "Data"))
K_modes_target <- as.integer(getOption("sst_spectral_ens_K_modes", getOption("sst_spectral_K_modes", 500L)))

problem_data <- build_sst_problem_data_ens(
  data_dir = data_dir,
  K_modes_target = K_modes_target
)

static <- problem_data$static
setup_snapshot <- problem_data$setup_snapshot
meta <- problem_data$meta

K_modes <- as.integer(meta$K_modes)

##########
# 3) MCMC configuration and initialization
##########

# Modes in calK_0 are allowed to toggle between beta (B_k=0) and X (B_k=1).
K0_target <- as.integer(getOption("sst_spectral_ens_K0", getOption("sst_spectral_K0", 150L)))
K0 <- max(1L, min(K0_target, K_modes - 1L))
calK_0 <- 2:(K0 + 1L)

p_k_const <- as.numeric(getOption("sst_spectral_ens_p_k", getOption("sst_spectral_p_k", 0.2)))
if (length(p_k_const) != 1L || !is.finite(p_k_const) || p_k_const <= 0 || p_k_const >= 1) {
  stop("sst_spectral_ens_p_k must be a single finite value strictly between 0 and 1.")
}
p_k <- rep(p_k_const, length(calK_0))
names(p_k) <- calK_0

B_k <- rbinom(n = length(calK_0), size = 1, prob = p_k)
names(B_k) <- calK_0

# Initial full allocation: mode 1 fixed in X; modes outside calK_0 fixed in X.
B_init <- rep(1L, K_modes)
names(B_init) <- seq_len(K_modes)
B_init[calK_0] <- B_k
B_init[1] <- 1L

delta1 <- as.numeric(getOption("sst_spectral_ens_delta1", getOption("sst_spectral_delta1", 0.5)))
delta01 <- as.numeric(getOption("sst_spectral_ens_delta01", getOption("sst_spectral_delta01", 1e-5)))
delta2 <- as.numeric(getOption("sst_spectral_ens_delta2", getOption("sst_spectral_delta2", 0.5)))
delta02 <- as.numeric(getOption("sst_spectral_ens_delta02", getOption("sst_spectral_delta02", 1e-5)))

sig2c_init <- as.numeric(getOption("sst_spectral_ens_sig2c_init", getOption("sst_spectral_sig2c_init", 0.2)))
sig2o_init <- as.numeric(getOption("sst_spectral_ens_sig2o_init", getOption("sst_spectral_sig2o_init", 0.2)))
if (!is.finite(sig2c_init) || sig2c_init <= 0) stop("sst_spectral_ens_sig2c_init must be > 0.")
if (!is.finite(sig2o_init) || sig2o_init <= 0) stop("sst_spectral_ens_sig2o_init must be > 0.")

# Inverse-gamma hyperparameters for variance updates.
gamma_o <- as.numeric(getOption("sst_spectral_ens_gamma_o", getOption("sst_spectral_gamma_o", 4.0)))
xi_o <- as.numeric(getOption("sst_spectral_ens_xi_o", getOption("sst_spectral_xi_o", 0.1)))
gamma_c <- as.numeric(getOption("sst_spectral_ens_gamma_c", getOption("sst_spectral_gamma_c", 4.0)))
xi_c <- as.numeric(getOption("sst_spectral_ens_xi_c", getOption("sst_spectral_xi_c", 0.1)))

##########
# 4) Runtime options and parameter/static objects
##########

n_iter <- max(1L, as.integer(getOption("sst_spectral_ens_n_iter", getOption("sst_spectral_n_iter", 300L))))
burn_in <- as.integer(getOption("sst_spectral_ens_burn_in", getOption("sst_spectral_burn_in", floor(n_iter / 4))))
burn_in <- max(0L, min(burn_in, n_iter - 1L))

progress_every <- max(1L, as.integer(getOption("sst_spectral_ens_progress_every", getOption("sst_spectral_progress_every", 20L))))
MCstore <- max(1L, as.integer(getOption("sst_spectral_ens_MCstore", getOption("sst_spectral_MCstore", progress_every))))
update_sigma2yo <- isTRUE(getOption(
  "sst_spectral_ens_update_sigma2yo",
  TRUE
))
update_sigma2yc <- isTRUE(getOption(
  "sst_spectral_ens_update_sigma2yc",
  TRUE
))
update_sigma2_after <- as.integer(getOption(
  "sst_spectral_ens_update_sigma2_after",
  15000L
))
if (length(update_sigma2_after) != 1L || is.na(update_sigma2_after)) {
  stop("sst_spectral_ens_update_sigma2_after must be a single integer >= 0.")
}
update_sigma2_after <- max(0L, update_sigma2_after)

rng_seed <- as.integer(getOption("sst_spectral_ens_seed", getOption("sst_spectral_seed", 123L)))
set.seed(rng_seed)

results_dir <- getOption(
  "sst_spectral_ens_results_dir",
  getOption("sst_spectral_results_dir", file.path("Data", "results"))
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

mcmc_rdata_file <- getOption(
  "sst_spectral_ens_mcmc_rdata",
  file.path(results_dir, sprintf("spectral2_sst_ens_mcmc_%s.RData", format(Sys.time(), "%Y%m%d_%H%M%S")))
)

params <- list(
  delta1 = delta1,
  delta01 = delta01,
  delta2 = delta2,
  delta02 = delta02,
  gamma_o = gamma_o,
  xi_o = xi_o,
  gamma_c = gamma_c,
  xi_c = xi_c,
  calK_0 = calK_0,
  p_k = p_k
)

setup_snapshot$seed <- rng_seed
setup_snapshot$K_modes <- K_modes
setup_snapshot$K0 <- K0
setup_snapshot$p_k_const <- p_k_const

##########
# 5) Run MCMC
##########

cat("Starting SST spectral ensemble MCMC with settings:\n")
cat(sprintf("  ocean sites (S): %d\n", static$S))
cat(sprintf("  observations (N): %d\n", static$N))
cat(sprintf("  ensemble size (M): %d\n", static$M_ens))
cat(sprintf("  within-ensemble sumsq (Wc): %.6f\n", static$Wc))
cat(sprintf("  K_modes: %d\n", K_modes))
cat(sprintf("  K0 (toggle set size): %d\n", length(calK_0)))
cat(sprintf("  delta1: %.6g, delta01: %.6g\n", delta1, delta01))
cat(sprintf("  delta2: %.6g, delta02: %.6g\n", delta2, delta02))
cat(sprintf("  n_iter: %d, burn_in: %d, MCstore: %d\n", n_iter, burn_in, MCstore))
cat(sprintf("  update_sigma2yo: %s\n", update_sigma2yo))
cat(sprintf("  update_sigma2yc: %s\n", update_sigma2yc))
cat(sprintf("  update_sigma2_after: %d\n", update_sigma2_after))
cat(sprintf("  output: %s\n", normalizePath(mcmc_rdata_file, mustWork = FALSE)))

mcmc_out <- run_mcmc_ens(
  B_init = B_init,
  sig2o_init = sig2o_init,
  sig2c_init = sig2c_init,
  n_iter = n_iter,
  burn_in = burn_in,
  progress_every = progress_every,
  MCstore = MCstore,
  static = static,
  params = params,
  update_sigma2yo = update_sigma2yo,
  update_sigma2yc = update_sigma2yc,
  update_sigma2_after = update_sigma2_after,
  save_rdata_path = mcmc_rdata_file,
  setup_snapshot = setup_snapshot
)

cat("MCMC finished.\n")
cat(sprintf("Saved artifact: %s\n", mcmc_out$save_rdata_path %||% "not saved"))
