update_B_local <- function(B_vec, sig2c_curr, sig2o_curr, static, params) {
  inv_sig2c <- 1 / sig2c_curr
  inv_sig2o <- 1 / sig2o_curr
  inv_sig2c2 <- inv_sig2c^2
  eps_num <- 1e-12

  I_X_work <- as.integer(which(B_vec == 1L))
  A_X_work <- static$U_obs[, I_X_work, drop = FALSE]
  D_X_work <- params$delta2 * static$lambda[I_X_work] + params$delta02
  Z_X_work <- static$Zc_all[I_X_work]
  h_X_work <- inv_sig2c * Z_X_work + inv_sig2o * as.numeric(crossprod(A_X_work, static$Yo_obs))

  P_X_work <- crossprod(A_X_work) * inv_sig2o
  if (length(I_X_work) > 0) {
    diag(P_X_work) <- diag(P_X_work) + D_X_work + inv_sig2c
  }

  for (k in params$calK_0) {
    pk_curr <- as.numeric(params$p_k[as.character(k)])
    pk_curr <- min(max(pk_curr, eps_num), 1 - eps_num)

    z_k <- static$Zc_all[k]
    d_X_k <- params$delta2 * static$lambda[k] + params$delta02
    d_beta_k <- params$delta1 * static$lambda[k] + params$delta01

    beta_term_k <- 0.5 * (
      log(d_beta_k) - log(d_beta_k + inv_sig2c) +
      (inv_sig2c2 * z_k^2) / (d_beta_k + inv_sig2c)
    )

    a_k_vec <- static$U_obs[, k]
    a_k_sq <- sum(a_k_vec^2)
    y_a_k <- sum(a_k_vec * static$Yo_obs)

    if (B_vec[k] == 1L) {
      idx_k <- match(k, I_X_work)
      a_kk <- P_X_work[idx_k, idx_k]
      h_k <- h_X_work[idx_k]

      if (length(I_X_work) > 1) {
        idx_r <- seq_along(I_X_work)[-idx_k]
        b_vec <- as.numeric(P_X_work[idx_r, idx_k, drop = FALSE])
        C_mat <- P_X_work[idx_r, idx_r, drop = FALSE]
        h_r <- h_X_work[idx_r]

        C_chol <- chol(C_mat)
        C_inv_b <- as.numeric(backsolve(C_chol, forwardsolve(t(C_chol), b_vec)))
        C_inv_hr <- as.numeric(backsolve(C_chol, forwardsolve(t(C_chol), h_r)))

        s_k <- as.numeric(a_kk - sum(b_vec * C_inv_b))
        g_k <- as.numeric(h_k - sum(b_vec * C_inv_hr))
      } else {
        s_k <- as.numeric(a_kk)
        g_k <- as.numeric(h_k)
      }

      s_k <- max(s_k, eps_num)
      delta_X_10 <- 0.5 * (-log(d_X_k) + log(s_k) - (g_k^2) / s_k)
      delta_Psi_10 <- log((1 - pk_curr) / pk_curr) + beta_term_k + delta_X_10

      prob_one <- inv_logit(-delta_Psi_10)
      B_new <- as.integer(rbinom(1, size = 1, prob = prob_one))

      if (B_new == 0L) {
        B_vec[k] <- 0L

        keep_idx <- seq_along(I_X_work)[-idx_k]
        I_X_work <- I_X_work[keep_idx]
        A_X_work <- A_X_work[, keep_idx, drop = FALSE]
        D_X_work <- D_X_work[keep_idx]
        Z_X_work <- Z_X_work[keep_idx]
        h_X_work <- h_X_work[keep_idx]

        if (length(keep_idx) > 0) {
          P_X_work <- P_X_work[keep_idx, keep_idx, drop = FALSE]
        } else {
          P_X_work <- matrix(0, nrow = 0, ncol = 0)
        }
      }
    } else {
      d_X_old <- length(I_X_work)
      a_add <- d_X_k + inv_sig2c + inv_sig2o * a_k_sq
      h_add <- inv_sig2c * z_k + inv_sig2o * y_a_k

      if (d_X_old > 0) {
        b_add <- inv_sig2o * as.numeric(crossprod(A_X_work, a_k_vec))
        C_chol <- chol(P_X_work)
        C_inv_b <- as.numeric(backsolve(C_chol, forwardsolve(t(C_chol), b_add)))
        C_inv_hr <- as.numeric(backsolve(C_chol, forwardsolve(t(C_chol), h_X_work)))

        s_k <- as.numeric(a_add - sum(b_add * C_inv_b))
        g_k <- as.numeric(h_add - sum(b_add * C_inv_hr))
      } else {
        b_add <- numeric(0)
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
        A_X_work <- cbind(A_X_work, a_k_vec)
        D_X_work <- c(D_X_work, d_X_k)
        Z_X_work <- c(Z_X_work, z_k)
        h_X_work <- c(h_X_work, h_add)

        if (d_X_old > 0) {
          P_X_work <- rbind(
            cbind(P_X_work, b_add),
            c(b_add, a_add)
          )
        } else {
          P_X_work <- matrix(a_add, nrow = 1, ncol = 1)
        }
      }
    }
  }

  B_vec
}

sample_b <- function(state, sig2c_curr) {
  d_beta <- state$d_beta

  if (d_beta > 0) {
    Pb_diag <- state$D_beta + (1 / sig2c_curr)
    rhs_b <- state$Zc_beta / sig2c_curr
    m_b <- rhs_b / Pb_diag
    b <- as.numeric(m_b + rnorm(d_beta) / sqrt(Pb_diag))
  } else {
    Pb_diag <- numeric(0)
    m_b <- numeric(0)
    b <- numeric(0)
  }

  list(b = b, m_b = m_b, Pb_diag = Pb_diag)
}

sample_x <- function(state, sig2c_curr, sig2o_curr, static) {
  d_X <- state$d_X
  A_X <- static$calH %*% state$U_X

  if (d_X > 0) {
    Px <- crossprod(A_X) / sig2o_curr
    diag(Px) <- diag(Px) + state$D_X + (1 / sig2c_curr)

    rhs_x <- state$Zc_X / sig2c_curr + as.numeric(crossprod(A_X, static$Yo)) / sig2o_curr
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

update_sigmas <- function(state, b, x, A_X, static, params) {
  X_mean_obs <- as.numeric(A_X %*% x)
  sse_o <- sum((static$Yo_obs - X_mean_obs)^2)
  sig2o_new <- rinvgamma1(params$gamma_o + static$N / 2, params$xi_o + sse_o / 2)

  beta_vec <- if (state$d_beta > 0) as.numeric(state$U_beta %*% b) else rep(0, static$S)
  X_vec <- as.numeric(state$U_X %*% x)
  sse_c <- sum((static$Yc - beta_vec - X_vec)^2)
  sig2c_new <- rinvgamma1(params$gamma_c + static$S / 2, params$xi_c + sse_c / 2)

  list(sig2o = sig2o_new, sig2c = sig2c_new, beta_vec = beta_vec, X_vec = X_vec)
}
