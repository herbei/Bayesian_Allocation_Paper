module_dir <- local({
  frame_files <- vapply(
    sys.frames(),
    function(env) {
      if (is.null(env$ofile)) {
        ""
      } else {
        env$ofile
      }
    },
    character(1L)
  )
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0L) {
    dirname(normalizePath(frame_files[length(frame_files)], winslash = "/", mustWork = FALSE))
  } else {
    cwd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    if (basename(cwd) == "modules") cwd else file.path(cwd, "modules")
  }
})
project_dir <- normalizePath(file.path(module_dir, ".."), winslash = "/", mustWork = FALSE)

source(file.path(module_dir, "mcmc_utils.R"))
source(file.path(module_dir, "mcmc_storage.R"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

read_grid_field <- function(path, n_rows, n_cols) {
  field_mat <- as.matrix(utils::read.table(path))
  if (!all(dim(field_mat) == c(n_rows, n_cols))) {
    stop(
      sprintf(
        "Expected %s to have dimensions %d x %d, found %d x %d.",
        path, n_rows, n_cols, nrow(field_mat), ncol(field_mat)
      ),
      call. = FALSE
    )
  }
  field_mat
}

build_rect_grid_laplacian <- function(n_rows, n_cols) {
  n_sites <- n_rows * n_cols
  L <- matrix(0, nrow = n_sites, ncol = n_sites)
  to_index <- function(r, c) r + (c - 1L) * n_rows

  for (c in seq_len(n_cols)) {
    for (r in seq_len(n_rows)) {
      s <- to_index(r, c)
      neighbors <- integer(0)
      if (r > 1L) {
        neighbors <- c(neighbors, to_index(r - 1L, c))
      }
      if (r < n_rows) {
        neighbors <- c(neighbors, to_index(r + 1L, c))
      }
      if (c > 1L) {
        neighbors <- c(neighbors, to_index(r, c - 1L))
      }
      if (c < n_cols) {
        neighbors <- c(neighbors, to_index(r, c + 1L))
      }
      L[s, s] <- length(neighbors)
      L[s, neighbors] <- -1
    }
  }

  L
}

load_or_compute_grid_eigendecomp <- function(eig_file, n_rows, n_cols) {
  if (file.exists(eig_file)) {
    eig_obj <- readRDS(eig_file)
    if (!is.null(eig_obj$U) && !is.null(eig_obj$eigenvalues)) {
      return(eig_obj)
    }
  }

  message(sprintf(
    "Computing full eigendecomposition for %d x %d OXY grid Laplacian.",
    n_rows, n_cols
  ))
  L <- build_rect_grid_laplacian(n_rows = n_rows, n_cols = n_cols)
  eig <- eigen(L, symmetric = TRUE)
  eig_obj <- list(
    U = eig$vectors,
    eigenvalues = eig$values,
    grid_dims = c(n_rows, n_cols)
  )
  saveRDS(eig_obj, eig_file)
  eig_obj
}

# Build OXY problem objects for the spectral MCMC workflow (ensemble version).
build_oxy_problem_data_ens <- function(data_dir = "Data", K_modes_target = 703L, M_ens_target = NULL,
                                       input_dir = file.path(data_dir, "inputs")) {
  n_rows <- 19L
  n_cols <- 37L
  n_sites <- n_rows * n_cols

  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = FALSE)
  legacy_input_files <- file.path(data_dir, c(
    "phi_o_data.txt",
    "api.txt",
    "oxy_ensemble_20_members.rds",
    "oxy_full_grid_laplacian_eigendecomp.rds"
  ))
  if (!dir.exists(input_dir) && all(file.exists(legacy_input_files))) {
    input_dir <- normalizePath(data_dir, winslash = "/", mustWork = FALSE)
  }

  yo_file <- file.path(input_dir, "phi_o_data.txt")
  api_file <- file.path(input_dir, "api.txt")
  ensemble_file <- file.path(input_dir, "oxy_ensemble_20_members.rds")
  eig_file <- file.path(input_dir, "oxy_full_grid_laplacian_eigendecomp.rds")

  required_files <- c(yo_file, api_file, ensemble_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop("Missing required OXY inputs:\n", paste0("  - ", missing_files, collapse = "\n"), call. = FALSE)
  }

  yo_mat <- read_grid_field(yo_file, n_rows = n_rows, n_cols = n_cols)
  yo_full <- as.vector(yo_mat)

  obs_sites <- scan(api_file, what = integer(), quiet = TRUE)
  if (length(obs_sites) == 0L) {
    stop("No observation indices found in api.txt.", call. = FALSE)
  }
  if (any(!is.finite(obs_sites)) || any(obs_sites < 1L) || any(obs_sites > n_sites)) {
    stop("Observation indices in api.txt must be integers between 1 and 703.", call. = FALSE)
  }

  ens <- readRDS(ensemble_file)
  required_ens_fields <- c("domain", "members", "fields")
  if (!all(required_ens_fields %in% names(ens))) {
    stop(
      "Ensemble RDS is missing required fields. Expected: ",
      paste(required_ens_fields, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(dim(ens$fields)) != 3L) {
    stop("Expected ensemble fields to be a 3D array [R x C x M].", call. = FALSE)
  }
  if (!identical(dim(ens$fields)[1:2], c(n_rows, n_cols))) {
    stop("OXY ensemble field dimensions do not match the expected 19 x 37 grid.", call. = FALSE)
  }

  member_labels_all <- as.character(ens$members)
  M_ens_available <- dim(ens$fields)[3]
  if (length(member_labels_all) != M_ens_available) {
    member_labels_all <- sprintf("m%02d", seq_len(M_ens_available))
  }

  if (is.null(M_ens_target) || length(M_ens_target) == 0L || is.na(M_ens_target[1L])) {
    M_ens_target <- M_ens_available
  }
  M_ens_target <- as.integer(M_ens_target[1L])
  if (!is.finite(M_ens_target) || M_ens_target < 1L || M_ens_target > M_ens_available) {
    stop(
      sprintf("Requested ensemble size M=%s must be an integer between 1 and %d.", M_ens_target, M_ens_available),
      call. = FALSE
    )
  }

  selected_member_indices <- if (M_ens_target < M_ens_available) {
    sort(sample.int(M_ens_available, size = M_ens_target, replace = FALSE))
  } else {
    seq_len(M_ens_available)
  }
  member_labels <- member_labels_all[selected_member_indices]
  M_ens <- length(selected_member_indices)
  ensemble_subsampled <- M_ens < M_ens_available

  Yc_full_matrix <- matrix(NA_real_, nrow = n_sites, ncol = M_ens)
  colnames(Yc_full_matrix) <- member_labels
  for (j in seq_len(M_ens)) {
    Yc_full_matrix[, j] <- as.vector(ens$fields[, , selected_member_indices[j]])
  }

  imputed_nonfinite_per_member <- integer(M_ens)
  imputed_fill_per_member <- rep(NA_real_, M_ens)
  for (j in seq_len(M_ens)) {
    yj <- as.numeric(Yc_full_matrix[, j])
    bad_j <- which(!is.finite(yj))
    imputed_nonfinite_per_member[j] <- length(bad_j)
    if (length(bad_j) > 0L) {
      fill_j <- stats::median(yj[is.finite(yj)])
      if (!is.finite(fill_j)) {
        fill_j <- 0
        warning(sprintf(
          "Member %s had no finite Yc values; imputed all with 0.",
          member_labels[j]
        ))
      }
      yj[bad_j] <- fill_j
      imputed_fill_per_member[j] <- fill_j
      Yc_full_matrix[, j] <- yj
    }
  }

  Yc_bar_full <- rowMeans(Yc_full_matrix)
  Wc <- sum((Yc_full_matrix - Yc_bar_full)^2)

  Yo_obs <- as.numeric(yo_full[obs_sites])
  if (any(!is.finite(Yo_obs))) {
    stop("Yo contains non-finite values at one or more api locations.", call. = FALSE)
  }
  N_obs <- length(Yo_obs)

  calH <- matrix(0, nrow = N_obs, ncol = n_sites)
  calH[cbind(seq_len(N_obs), obs_sites)] <- 1

  eig_obj <- load_or_compute_grid_eigendecomp(
    eig_file = eig_file,
    n_rows = n_rows,
    n_cols = n_cols
  )
  if (is.null(eig_obj$U) || is.null(eig_obj$eigenvalues)) {
    stop("Eigendecomposition object must contain U and eigenvalues.", call. = FALSE)
  }

  lambda_all <- as.numeric(eig_obj$eigenvalues)
  if (length(lambda_all) != n_sites) {
    stop("Eigenvalue count does not match the number of OXY grid sites.", call. = FALSE)
  }
  if (!all(dim(eig_obj$U) == c(n_sites, n_sites))) {
    stop("Eigenvector matrix dimension mismatch for OXY eigendecomposition.", call. = FALSE)
  }

  ord <- order(lambda_all, decreasing = FALSE)
  K_modes <- max(2L, min(as.integer(K_modes_target), n_sites))
  mode_indices <- ord[seq_len(K_modes)]

  lambda_modes <- lambda_all[mode_indices]
  U_modes <- eig_obj$U[, mode_indices, drop = FALSE]

  U_obs <- U_modes[obs_sites, , drop = FALSE]
  G_obs <- crossprod(U_obs)
  Yo_proj_all <- as.numeric(crossprod(U_obs, Yo_obs))

  static <- list(
    S = n_sites,
    R = n_rows,
    C = n_cols,
    U = U_modes,
    lambda = lambda_modes,
    Yc_bar = as.numeric(Yc_bar_full),
    M_ens = as.integer(M_ens),
    Wc = as.numeric(Wc),
    Yo_obs = as.numeric(Yo_obs),
    calH = calH,
    N = N_obs,
    U_obs = U_obs,
    G_obs = G_obs,
    U_obs_col_sq = as.numeric(diag(G_obs)),
    Yo_proj_all = Yo_proj_all,
    Zc_all = as.numeric(crossprod(U_modes, Yc_bar_full))
  )

  setup_snapshot <- list(
    data_dir = normalizePath(data_dir, mustWork = FALSE),
    input_dir = normalizePath(input_dir, mustWork = FALSE),
    grid_dims = list(R = n_rows, C = n_cols, S_full = n_sites),
    obs_count = N_obs,
    obs_sites = obs_sites,
    mode_indices_in_full_basis = mode_indices,
    ensemble_count = M_ens,
    ensemble_count_requested = M_ens_target,
    ensemble_count_available = M_ens_available,
    ensemble_subsampled = ensemble_subsampled,
    ensemble_member_indices_selected = selected_member_indices,
    ensemble_member_labels_selected = member_labels,
    ensemble_member_labels = member_labels,
    Wc = as.numeric(Wc),
    imputed_nonfinite_per_member = imputed_nonfinite_per_member,
    imputed_fill_per_member = imputed_fill_per_member,
    eig_file = normalizePath(eig_file, mustWork = FALSE)
  )

  meta <- list(
    K_modes = K_modes,
    N_obs = N_obs,
    R = n_rows,
    C = n_cols,
    S_full = n_sites,
    M_ens = M_ens,
    M_ens_available = M_ens_available
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

      if (d_X_old > 0L) {
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

        if (d_X_old > 0L) {
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
    stop(sprintf("Length of B_vec (%d) must match number of modes (%d).", length(B_vec), n_modes), call. = FALSE)
  }

  I_beta_loc <- as.integer(which(B_vec == 0L))
  I_X_loc <- as.integer(setdiff(seq_len(n_modes), I_beta_loc))

  d_beta_loc <- length(I_beta_loc)
  d_X_loc <- length(I_X_loc)

  U_beta_loc <- static$U[, I_beta_loc, drop = FALSE]
  U_X_loc <- static$U[, I_X_loc, drop = FALSE]

  D_beta_loc <- params$delta1 * static$lambda[I_beta_loc] + params$delta01
  D_X_loc <- params$delta2 * static$lambda[I_X_loc] + params$delta02

  Zc_beta_loc <- if (d_beta_loc > 0L) as.numeric(static$Zc_all[I_beta_loc]) else numeric(0)
  Zc_X_loc <- if (d_X_loc > 0L) as.numeric(static$Zc_all[I_X_loc]) else numeric(0)

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

  if (d_beta > 0L) {
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
  A_X <- static$U_obs[, state$I_X, drop = FALSE]
  inv_sig2c_eff <- static$M_ens / sig2c_curr

  if (d_X > 0L) {
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

  beta_vec <- if (state$d_beta > 0L) as.numeric(state$U_beta %*% b) else rep(0, static$S)
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

  if (state$d_X > 0L) {
    resid_ZXc <- state$Zc_X - x_curr
    logp_ZXc <- 0.5 * state$d_X * log(two_pi * sig2c_eff) + 0.5 * sum(resid_ZXc^2) / sig2c_eff
  } else {
    logp_ZXc <- 0
  }

  if (state$d_beta > 0L) {
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
    stop("update_sigma2_after must be a single non-missing integer >= 0.", call. = FALSE)
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
    B_curr <- update_B_local_ens(B_curr, sig2c_curr, sig2o_curr, static, params)
    state_curr <- build_state_from_B_ens(B_curr, static, params)

    b_res <- sample_b_ens(state_curr, sig2c_curr, static)
    b_curr <- b_res$b

    x_res <- sample_x_ens(state_curr, sig2c_curr, sig2o_curr, static)
    x_curr <- x_res$x

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
      beta_vec_curr <- if (state_curr$d_beta > 0L) as.numeric(state_curr$U_beta %*% b_curr) else rep(0, static$S)
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

    if (iter %% progress_every == 0L || iter == 1L || iter == n_iter) {
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
# 2) Data loading and graph-basis setup (OXY full grid, ensemble Yc)
##########

data_dir <- getOption("oxy_spectral_ens_data_dir", file.path(project_dir, "Data"))
input_dir <- getOption("oxy_spectral_ens_input_dir", file.path(data_dir, "inputs"))
K_modes_target <- as.integer(getOption("oxy_spectral_ens_K_modes", 703L))
rng_seed <- as.integer(getOption("oxy_spectral_ens_seed", 123L))
M_ens_target <- as.integer(getOption("oxy_spectral_ens_M_ens", getOption("oxy_spectral_ens_M", 20L)))

if (length(rng_seed) != 1L || is.na(rng_seed)) {
  stop("oxy_spectral_ens_seed must be a single integer.", call. = FALSE)
}
if (length(M_ens_target) != 1L || is.na(M_ens_target)) {
  stop("oxy_spectral_ens_M_ens must be a single integer.", call. = FALSE)
}

set.seed(rng_seed)

problem_data <- build_oxy_problem_data_ens(
  data_dir = data_dir,
  input_dir = input_dir,
  K_modes_target = K_modes_target,
  M_ens_target = M_ens_target
)

static <- problem_data$static
setup_snapshot <- problem_data$setup_snapshot
meta <- problem_data$meta

K_modes <- as.integer(meta$K_modes)

##########
# 3) MCMC configuration and initialization
##########

K0_target <- as.integer(getOption("oxy_spectral_ens_K0", 300L))
K0 <- max(1L, min(K0_target, K_modes - 1L))
calK_0 <- 2:(K0 + 1L)

p_k_const <- as.numeric(getOption("oxy_spectral_ens_p_k", 0.5))
if (length(p_k_const) != 1L || !is.finite(p_k_const) || p_k_const <= 0 || p_k_const >= 1) {
  stop("oxy_spectral_ens_p_k must be a single finite value strictly between 0 and 1.", call. = FALSE)
}
p_k <- rep(p_k_const, length(calK_0))
names(p_k) <- calK_0

B_k <- rbinom(n = length(calK_0), size = 1, prob = p_k)
names(B_k) <- calK_0

B_init <- rep(1L, K_modes)
names(B_init) <- seq_len(K_modes)
B_init[calK_0] <- B_k
B_init[1] <- 1L
B_init_source <- "random prior draw"

B_init_file <- getOption("oxy_spectral_ens_B_init_file", file.path(input_dir, "Binit.txt"))
if (length(B_init_file) != 1L || is.na(B_init_file)) {
  stop("oxy_spectral_ens_B_init_file must be a single path or an empty string.", call. = FALSE)
}
if (nzchar(B_init_file)) {
  if (!file.exists(B_init_file)) {
    stop(sprintf("B_init file does not exist: %s", B_init_file), call. = FALSE)
  }
  B_init_from_file <- scan(B_init_file, what = numeric(), quiet = TRUE)
  if (length(B_init_from_file) != K_modes) {
    stop(
      sprintf(
        "B_init file %s has length %d; expected %d.",
        B_init_file, length(B_init_from_file), K_modes
      ),
      call. = FALSE
    )
  }
  if (any(!is.finite(B_init_from_file)) || any(!(B_init_from_file %in% c(0, 1)))) {
    stop("B_init file must contain only binary 0/1 values.", call. = FALSE)
  }

  B_init <- as.integer(B_init_from_file)
  names(B_init) <- seq_len(K_modes)
  fixed_one_modes <- setdiff(seq_len(K_modes), calK_0)
  if (any(B_init[fixed_one_modes] != 1L)) {
    warning("B_init file has zeros outside calK_0; fixed non-toggle modes were reset to 1.")
    B_init[fixed_one_modes] <- 1L
  }
  B_init[1] <- 1L
  B_init_source <- normalizePath(B_init_file, mustWork = FALSE)
}

delta_common <- as.numeric(getOption("oxy_spectral_ens_delta", 0.5))
delta1 <- as.numeric(getOption("oxy_spectral_ens_delta1", delta_common))
delta01 <- as.numeric(getOption("oxy_spectral_ens_delta01", 1e-5))
delta2 <- as.numeric(getOption("oxy_spectral_ens_delta2", delta_common))
delta02 <- as.numeric(getOption("oxy_spectral_ens_delta02", 1e-5))
if (length(delta1) != 1L || !is.finite(delta1) || delta1 <= 0) {
  stop("oxy_spectral_ens_delta1 must be a single finite positive value.", call. = FALSE)
}
if (length(delta01) != 1L || !is.finite(delta01) || delta01 <= 0) {
  stop("oxy_spectral_ens_delta01 must be a single finite positive value.", call. = FALSE)
}
if (length(delta2) != 1L || !is.finite(delta2) || delta2 <= 0) {
  stop("oxy_spectral_ens_delta2 must be a single finite positive value.", call. = FALSE)
}
if (length(delta02) != 1L || !is.finite(delta02) || delta02 <= 0) {
  stop("oxy_spectral_ens_delta02 must be a single finite positive value.", call. = FALSE)
}

sig2c_init <- as.numeric(getOption("oxy_spectral_ens_sig2c_init", 0.2))
sig2o_init <- as.numeric(getOption("oxy_spectral_ens_sig2o_init", 0.2))
if (!is.finite(sig2c_init) || sig2c_init <= 0) stop("oxy_spectral_ens_sig2c_init must be > 0.", call. = FALSE)
if (!is.finite(sig2o_init) || sig2o_init <= 0) stop("oxy_spectral_ens_sig2o_init must be > 0.", call. = FALSE)

gamma_o <- as.numeric(getOption("oxy_spectral_ens_gamma_o", 4.0))
xi_o <- as.numeric(getOption("oxy_spectral_ens_xi_o", 0.1))
gamma_c <- as.numeric(getOption("oxy_spectral_ens_gamma_c", 4.0))
xi_c <- as.numeric(getOption("oxy_spectral_ens_xi_c", 0.1))

##########
# 4) Runtime options and parameter/static objects
##########

n_iter <- max(1L, as.integer(getOption("oxy_spectral_ens_n_iter", 300L)))
burn_in <- as.integer(getOption("oxy_spectral_ens_burn_in", floor(n_iter / 4)))
burn_in <- max(0L, min(burn_in, n_iter - 1L))

progress_every <- max(1L, as.integer(getOption("oxy_spectral_ens_progress_every", 20L)))
MCstore <- max(1L, as.integer(getOption("oxy_spectral_ens_MCstore", progress_every)))
update_sigma2yo <- isTRUE(getOption("oxy_spectral_ens_update_sigma2yo", TRUE))
update_sigma2yc <- isTRUE(getOption("oxy_spectral_ens_update_sigma2yc", TRUE))
update_sigma2_after <- as.integer(getOption("oxy_spectral_ens_update_sigma2_after", 15000L))
if (length(update_sigma2_after) != 1L || is.na(update_sigma2_after)) {
  stop("oxy_spectral_ens_update_sigma2_after must be a single integer >= 0.", call. = FALSE)
}
update_sigma2_after <- max(0L, update_sigma2_after)

set.seed(rng_seed)

results_dir <- file.path(data_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

mcmc_rdata_file <- getOption(
  "oxy_spectral_ens_mcmc_rdata",
  file.path(results_dir, sprintf("spectral2_ens_mcmc_%s.RData", format(Sys.time(), "%Y%m%d_%H%M%S")))
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
setup_snapshot$B_init_source <- B_init_source

##########
# 5) Run MCMC
##########

cat("Starting OXY spectral ensemble MCMC with settings:\n")
cat(sprintf("  grid sites (S): %d\n", static$S))
cat(sprintf("  observations (N): %d\n", static$N))
cat(sprintf("  ensemble size (M): %d of %d available\n", static$M_ens, setup_snapshot$ensemble_count_available %||% static$M_ens))
if (isTRUE(setup_snapshot$ensemble_subsampled)) {
  cat(sprintf("  selected ensemble members: %s\n", paste(setup_snapshot$ensemble_member_labels_selected, collapse = ", ")))
}
cat(sprintf("  within-ensemble sumsq (Wc): %.6f\n", static$Wc))
cat(sprintf("  K_modes: %d\n", K_modes))
cat(sprintf("  K0 (toggle set size): %d\n", length(calK_0)))
cat(sprintf("  B_init source: %s\n", B_init_source))
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
