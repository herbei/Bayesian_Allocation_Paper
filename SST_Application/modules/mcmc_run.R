compute_log_mcmc_summaries <- function(state, x_curr, b_curr, A_X, static, sig2o_curr, sig2c_curr) {
  two_pi <- 2 * pi

  # -log p(Yo | x, B, sigma^2_yo) with Yo observed at calH rows.
  X_mean_obs <- as.numeric(A_X %*% x_curr)
  resid_o <- static$Yo_obs - X_mean_obs
  logp_Yo <- 0.5 * static$N * log(two_pi * sig2o_curr) + 0.5 * sum(resid_o^2) / sig2o_curr

  # -log p(Z_X^c | x, B, sigma^2_yc) in the transformed X block.
  if (state$d_X > 0) {
    resid_ZXc <- state$Zc_X - x_curr
    logp_ZXc <- 0.5 * state$d_X * log(two_pi * sig2c_curr) + 0.5 * sum(resid_ZXc^2) / sig2c_curr
  } else {
    logp_ZXc <- 0
  }

  # -log p(Z_beta^c | b, B, sigma^2_yo) as requested.
  if (state$d_beta > 0) {
    resid_Zbetac <- state$Zc_beta - b_curr
    logp_Zbetac <- 0.5 * state$d_beta * log(two_pi * sig2o_curr) + 0.5 * sum(resid_Zbetac^2) / sig2o_curr
  } else {
    logp_Zbetac <- 0
  }

  list(
    logp_Yo = logp_Yo,
    logp_ZXc = logp_ZXc,
    logp_Zbetac = logp_Zbetac
  )
}

run_mcmc <- function(B_init, sig2o_init, sig2c_init, n_iter, burn_in,
                     progress_every, MCstore, static, params,
                     update_sigma2yo = TRUE,
                     update_sigma2yc = TRUE,
                     save_rdata_path = NULL,
                     setup_snapshot = NULL) {
  store <- init_storage(n_iter, MCstore, static, params)

  B_curr <- B_init
  sig2o_curr <- sig2o_init
  sig2c_curr <- sig2c_init

  state_curr <- build_state_from_B(B_curr, static, params)
  b_curr <- rep(0, state_curr$d_beta)
  x_curr <- rep(0, state_curr$d_X)

  for (iter in seq_len(n_iter)) {
    # Step 1: local one-mode collapsed update for B
    B_curr <- update_B_local(B_curr, sig2c_curr, sig2o_curr, static, params)

    # Step 2: rebuild state under updated B
    state_curr <- build_state_from_B(B_curr, static, params)

    # Step 3: sample b | B, sig2c, Yc
    b_res <- sample_b(state_curr, sig2c_curr)
    b_curr <- b_res$b

    # Step 4: sample x | B, sig2c, sig2o, Yc, Yo
    x_res <- sample_x(state_curr, sig2c_curr, sig2o_curr, static)
    x_curr <- x_res$x

    # Step 5/6: sample sig2o and/or sig2c (or keep fixed if disabled)
    if (update_sigma2yo || update_sigma2yc) {
      sig_res <- update_sigmas(state_curr, b_curr, x_curr, x_res$A_X, static, params)
      if (update_sigma2yo) {
        sig2o_curr <- sig_res$sig2o
      }
      if (update_sigma2yc) {
        sig2c_curr <- sig_res$sig2c
      }
      beta_vec_curr <- sig_res$beta_vec
      X_vec_curr <- sig_res$X_vec
    } else {
      beta_vec_curr <- if (state_curr$d_beta > 0) as.numeric(state_curr$U_beta %*% b_curr) else rep(0, static$S)
      X_vec_curr <- as.numeric(state_curr$U_X %*% x_curr)
    }

    log_summ <- compute_log_mcmc_summaries(
      state = state_curr,
      x_curr = x_curr,
      b_curr = b_curr,
      A_X = x_res$A_X,
      static = static,
      sig2o_curr = sig2o_curr,
      sig2c_curr = sig2c_curr
    )

    # Step 7: store if requested by MCstore
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
        source = "run_mcmc",
        format_version = 1L
      ),
      run_config = list(
        n_iter = n_iter,
        burn_in = burn_in,
        progress_every = progress_every,
        MCstore = MCstore,
        update_sigma2yo = update_sigma2yo,
        update_sigma2yc = update_sigma2yc,
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
