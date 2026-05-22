#!/usr/bin/env Rscript

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
})
setwd(script_dir)

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  as.integer(value)
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  as.numeric(value)
}

n_iter <- env_int("N_ITER", 30000L)
burn_in <- env_int("BURN_IN", floor(n_iter / 4))
progress_every <- env_int("PROGRESS_EVERY", 50L)
mcstore <- env_int("MCSTORE", 50L)
seed <- env_int("SEED", 1002L)

dir.create(file.path("Data", "results"), recursive = TRUE, showWarnings = FALSE)
output_file <- Sys.getenv(
  "OUTPUT_FILE",
  unset = file.path("Data", "results", "simulation_study_selected_run.RData")
)

options(
  oxy_spectral_ens_data_dir = "Data",
  oxy_spectral_ens_input_dir = file.path("Data", "inputs"),
  oxy_spectral_ens_K_modes = 703L,
  oxy_spectral_ens_K0 = env_int("K0", 300L),
  oxy_spectral_ens_M_ens = env_int("M_ENS", 20L),
  oxy_spectral_ens_p_k = env_num("P_K", 0.5),
  oxy_spectral_ens_delta1 = env_num("DELTA1", 2.5),
  oxy_spectral_ens_delta01 = env_num("DELTA01", 1e-4),
  oxy_spectral_ens_delta2 = env_num("DELTA2", 2.5),
  oxy_spectral_ens_delta02 = env_num("DELTA02", 1e-4),
  oxy_spectral_ens_sig2c_init = env_num("SIG2C_INIT", 0.5),
  oxy_spectral_ens_sig2o_init = env_num("SIG2O_INIT", 0.5),
  oxy_spectral_ens_gamma_o = 4.0,
  oxy_spectral_ens_xi_o = 0.1,
  oxy_spectral_ens_gamma_c = 4.0,
  oxy_spectral_ens_xi_c = 0.1,
  oxy_spectral_ens_n_iter = n_iter,
  oxy_spectral_ens_burn_in = burn_in,
  oxy_spectral_ens_progress_every = progress_every,
  oxy_spectral_ens_MCstore = mcstore,
  oxy_spectral_ens_update_sigma2yo = TRUE,
  oxy_spectral_ens_update_sigma2yc = TRUE,
  oxy_spectral_ens_update_sigma2_after = env_int("UPDATE_SIGMA2_AFTER", burn_in),
  oxy_spectral_ens_seed = seed,
  oxy_spectral_ens_mcmc_rdata = output_file
)

source(file.path("modules", "spectral2_ens.R"))
