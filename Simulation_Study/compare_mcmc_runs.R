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

mcmc_runs_dir <- "Data/results"

source(file.path("modules", "mcmc_run_comparison.R"))

compare_mcmc_runs(
  runs_dir = mcmc_runs_dir,
  invocation_dir = script_dir
)
