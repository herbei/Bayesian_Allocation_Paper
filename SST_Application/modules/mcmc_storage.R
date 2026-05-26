init_storage <- function(n_iter, MCstore, static, params) {
  if (MCstore > n_iter) {
    store_iters <- n_iter
  } else {
    store_iters <- seq.int(from = MCstore, to = n_iter, by = MCstore)
  }

  n_store <- length(store_iters)
  store_mask <- rep(FALSE, n_iter)
  store_mask[store_iters] <- TRUE

  list(
    n_store = n_store,
    store_iters = store_iters,
    store_mask = store_mask,
    store_pos = 0L,
    B_chain = matrix(NA_integer_, nrow = n_store, ncol = length(params$calK_0),
                     dimnames = list(NULL, paste0("k", params$calK_0))),
    sig2o_chain = numeric(n_store),
    sig2c_chain = numeric(n_store),
    eta1_chain = numeric(n_store),
    eta2_chain = numeric(n_store),
    delta1_chain = numeric(n_store),
    delta2_chain = numeric(n_store),
    logp_Yo_chain = numeric(n_store),
    logp_ZXc_chain = numeric(n_store),
    logp_Zbetac_chain = numeric(n_store),
    d_beta_chain = integer(n_store),
    iter_chain = integer(n_store),
    beta_chain = matrix(NA_real_, nrow = static$S, ncol = n_store),
    X_chain = matrix(NA_real_, nrow = static$S, ncol = n_store),
    b_chain = vector("list", n_store),
    x_chain = vector("list", n_store)
  )
}

store_draw <- function(store, iter, B_vec, state, sig2o_curr, sig2c_curr,
                       eta1_curr = NA_real_, eta2_curr = NA_real_,
                       delta1_curr = NA_real_, delta2_curr = NA_real_,
                       logp_Yo_curr, logp_ZXc_curr, logp_Zbetac_curr,
                       beta_vec, X_vec, b_vec, x_vec, params) {
  if (!store$store_mask[iter]) {
    return(store)
  }

  pos <- store$store_pos + 1L
  store$store_pos <- pos
  store$B_chain[pos, ] <- B_vec[params$calK_0]
  store$sig2o_chain[pos] <- sig2o_curr
  store$sig2c_chain[pos] <- sig2c_curr
  store$eta1_chain[pos] <- eta1_curr
  store$eta2_chain[pos] <- eta2_curr
  store$delta1_chain[pos] <- delta1_curr
  store$delta2_chain[pos] <- delta2_curr
  store$logp_Yo_chain[pos] <- logp_Yo_curr
  store$logp_ZXc_chain[pos] <- logp_ZXc_curr
  store$logp_Zbetac_chain[pos] <- logp_Zbetac_curr
  store$d_beta_chain[pos] <- state$d_beta
  store$iter_chain[pos] <- iter
  store$beta_chain[, pos] <- beta_vec
  store$X_chain[, pos] <- X_vec
  store$b_chain[[pos]] <- b_vec
  store$x_chain[[pos]] <- x_vec

  store
}
