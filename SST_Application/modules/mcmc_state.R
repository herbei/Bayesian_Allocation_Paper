build_state_from_B <- function(B_vec, static, params) {
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

  tilde_U_loc <- cbind(U_beta_loc, U_X_loc)
  Z_c_loc <- as.numeric(crossprod(tilde_U_loc, static$Yc))

  Zc_beta_loc <- if (d_beta_loc > 0) Z_c_loc[seq_len(d_beta_loc)] else numeric(0)
  Zc_X_loc <- if (d_X_loc > 0) Z_c_loc[(d_beta_loc + 1):(d_beta_loc + d_X_loc)] else numeric(0)

  list(
    I_beta = I_beta_loc,
    I_X = I_X_loc,
    d_beta = d_beta_loc,
    d_X = d_X_loc,
    U_beta = U_beta_loc,
    U_X = U_X_loc,
    D_beta = D_beta_loc,
    D_X = D_X_loc,
    tilde_U = tilde_U_loc,
    Z_c = Z_c_loc,
    Zc_beta = Zc_beta_loc,
    Zc_X = Zc_X_loc
  )
}
