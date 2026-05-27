`%||%` <- function(a, b) if (!is.null(a)) a else b

resolve_path <- function(path, base_dir) {
  path <- path.expand(path)
  if (grepl("^(/|[A-Za-z]:[\\\\/])", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path(base_dir, path), winslash = "/", mustWork = FALSE)
}

format_decimal <- function(x, digits = 3L) {
  x <- as.numeric(x)
  out <- rep("NA", length(x))
  ok <- is.finite(x)
  out[ok] <- formatC(x[ok], format = "f", digits = digits)
  out
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  gsub("\\^", "\\\\textasciicircum{}", x)
}

read_settings_file <- function(result_file) {
  settings_file <- sub("\\.RData$", ".sett.txt", result_file)
  if (!file.exists(settings_file)) return(list())

  lines <- readLines(settings_file, warn = FALSE)
  lines <- lines[nzchar(lines) & grepl("=", lines, fixed = TRUE)]
  if (length(lines) == 0L) return(list())

  keys <- sub("=.*$", "", lines)
  vals <- substring(lines, nchar(keys) + 2L)
  as.list(stats::setNames(vals, keys))
}

setting_num <- function(settings, key) {
  value <- settings[[key]]
  if (is.null(value) || length(value) == 0L || !nzchar(value[1L])) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(value[1L]))
}

setting_int <- function(settings, key) {
  value <- settings[[key]]
  if (is.null(value) || length(value) == 0L || !nzchar(value[1L])) {
    return(NA_integer_)
  }
  suppressWarnings(as.integer(value[1L]))
}

safe_rmse <- function(truth, estimate) {
  ok <- is.finite(truth) & is.finite(estimate)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((estimate[ok] - truth[ok])^2))
}

posterior_mean <- function(chain_matrix, keep_idx, fallback) {
  chain_matrix <- as.matrix(chain_matrix %||% matrix(numeric(0), nrow = 0L, ncol = 0L))
  keep_idx <- keep_idx[keep_idx <= ncol(chain_matrix)]
  if (length(keep_idx) == 0L) return(as.numeric(fallback %||% numeric(0)))
  rowMeans(chain_matrix[, keep_idx, drop = FALSE], na.rm = TRUE)
}

kept_mean <- function(x, keep_idx) {
  x <- as.numeric(x %||% numeric(0))
  keep_idx <- keep_idx[keep_idx <= length(x)]
  if (length(keep_idx) == 0L) return(NA_real_)
  mean(x[keep_idx], na.rm = TRUE)
}

run_id_from_file <- function(file_name) {
  match <- regmatches(file_name, regexpr("[0-9]{8}_[0-9]{6}_[0-9]{6}", file_name))
  if (length(match) == 0L || !nzchar(match)) sub("\\.RData$", "", file_name) else match
}

summarize_mcmc_file <- function(result_file, true_x) {
  result_env <- new.env(parent = emptyenv())
  load(result_file, envir = result_env)
  if (!exists("mcmc_artifact", envir = result_env, inherits = FALSE)) {
    stop("missing mcmc_artifact")
  }

  artifact <- get("mcmc_artifact", envir = result_env, inherits = FALSE)
  output <- artifact$output %||% list()
  store <- output$store %||% list()
  params <- artifact$params %||% list()
  setup <- artifact$setup %||% list()
  run_config <- artifact$run_config %||% list()
  static <- artifact$static %||% list()
  settings <- read_settings_file(result_file)

  Yc_bar <- as.numeric(static$Yc_bar %||% numeric(0))
  if (length(Yc_bar) == 0L) stop("missing static$Yc_bar")
  if (length(Yc_bar) != length(true_x)) {
    stop(sprintf("Yc_bar length %d does not match trueXX length %d", length(Yc_bar), length(true_x)))
  }

  iter_chain <- as.integer(store$iter_chain %||% integer(0))
  keep_idx <- as.integer(output$keep_idx %||% seq_along(iter_chain))
  keep_idx <- unique(keep_idx[is.finite(keep_idx) & keep_idx >= 1L])
  if (length(iter_chain) > 0L) keep_idx <- keep_idx[keep_idx <= length(iter_chain)]
  if (length(keep_idx) == 0L) stop("no kept MCMC draws available")

  burn_in <- as.integer(run_config$burn_in %||% 0L)
  update_sigma2_after <- as.integer(run_config$update_sigma2_after %||% burn_in)
  summary_after <- max(c(0L, burn_in, update_sigma2_after), na.rm = TRUE)
  comparison_keep_idx <- if (length(iter_chain) > 0L) keep_idx[iter_chain[keep_idx] > summary_after] else keep_idx
  if (length(comparison_keep_idx) == 0L) comparison_keep_idx <- tail(keep_idx, 1L)

  X_post_mean <- posterior_mean(store$X_chain, comparison_keep_idx, output$X_post_mean)
  beta_post_mean <- posterior_mean(store$beta_chain, comparison_keep_idx, output$beta_post_mean)
  if (length(X_post_mean) != length(true_x)) stop("X posterior mean length does not match trueXX")
  if (length(beta_post_mean) != length(true_x)) stop("beta posterior mean length does not match trueXX")

  file_name <- basename(result_file)
  run_id <- run_id_from_file(file_name)
  true_beta <- Yc_bar - true_x

  data.frame(
    file = file_name,
    run = run_id,
    timestamp = run_id,
    K0 = as.integer(setup$K0 %||% setting_int(settings, "k0")),
    M_ens = as.integer(setup$ensemble_count %||% setting_int(settings, "m_ens")),
    p_k = as.numeric(unique(params$p_k %||% setup$p_k_const %||% setting_num(settings, "p_k"))[1L]),
    seed = as.integer(setup$seed %||% setting_int(settings, "seed")),
    delta_1 = as.numeric(params$delta1 %||% setting_num(settings, "delta1")),
    delta_2 = as.numeric(params$delta2 %||% setting_num(settings, "delta2")),
    rmse_x = safe_rmse(true_x, X_post_mean),
    rmse_beta = safe_rmse(true_beta, beta_post_mean),
    post_mean_sigma2yc = kept_mean(store$sig2c_chain, comparison_keep_idx),
    post_mean_sigma2yo = kept_mean(store$sig2o_chain, comparison_keep_idx),
    post_mean_d_beta = kept_mean(store$d_beta_chain, comparison_keep_idx),
    n_kept = length(comparison_keep_idx),
    sig2yc_init = as.numeric(run_config$sig2c_init %||% setting_num(settings, "SIG2C")),
    sig2yo_init = as.numeric(run_config$sig2o_init %||% setting_num(settings, "SIG2O")),
    stringsAsFactors = FALSE
  )
}

write_comparison_longtable <- function(comparison_df, output_file, caption, label) {
  header_cells <- c(
    "Rank", "RMSE$(X)$", "RMSE$(\\beta)$", "Rank$(\\beta)$",
    "$\\overline{\\sigma^2_{yo}}$", "$\\overline{\\sigma^2_{yc}}$",
    "$p_k$", "$\\overline{d_{\\beta}}$", "$\\delta_1$", "$\\delta_2$"
  )

  row_lines <- vapply(seq_len(nrow(comparison_df)), function(i) {
    row <- comparison_df[i, ]
    cells <- c(
      as.character(row$rank_rmse_x),
      format_decimal(row$rmse_x, 3L),
      format_decimal(row$rmse_beta, 3L),
      as.character(row$rank_rmse_beta),
      format_decimal(row$post_mean_sigma2yo, 3L),
      format_decimal(row$post_mean_sigma2yc, 3L),
      format_decimal(row$p_k, 2L),
      format_decimal(row$post_mean_d_beta, 2L),
      format_decimal(row$delta_1, 2L),
      format_decimal(row$delta_2, 2L)
    )
    paste0(paste(cells, collapse = " & "), " \\\\")
  }, character(1L))

  writeLines(c(
    "% Requires \\usepackage{booktabs,longtable}",
    "% Generated by compare_mcmc_runs.R",
    "\\begin{longtable}{@{}rrrrrrrrrr@{}}",
    paste0("\\caption{", latex_escape(caption), "}\\label{", label, "}\\\\"),
    "\\toprule",
    paste0(paste(header_cells, collapse = " & "), " \\\\"),
    "\\midrule",
    "\\endfirsthead",
    "\\multicolumn{10}{l}{\\tablename\\ \\thetable\\ -- continued from previous page}\\\\",
    "\\toprule",
    paste0(paste(header_cells, collapse = " & "), " \\\\"),
    "\\midrule",
    "\\endhead",
    "\\midrule",
    "\\multicolumn{10}{r}{Continued on next page}\\\\",
    "\\midrule",
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot",
    row_lines,
    "\\end{longtable}"
  ), output_file)
}

compare_mcmc_runs <- function(runs_dir = "Data/results",
                              output_dir = NULL,
                              true_file = file.path("Data", "inputs", "trueXX.txt"),
                              invocation_dir = normalizePath(getwd(), winslash = "/", mustWork = FALSE)) {
  invocation_dir <- normalizePath(invocation_dir, winslash = "/", mustWork = FALSE)
  runs_dir <- resolve_path(runs_dir, invocation_dir)
  true_file <- resolve_path(true_file, invocation_dir)
  output_dir <- if (is.null(output_dir)) file.path(runs_dir, "comparison") else resolve_path(output_dir, invocation_dir)

  if (!dir.exists(runs_dir)) stop("MCMC runs directory does not exist: ", runs_dir, call. = FALSE)
  if (!file.exists(true_file)) stop("trueXX file does not exist: ", true_file, call. = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  result_files <- sort(normalizePath(
    list.files(runs_dir, pattern = "\\.RData$", full.names = TRUE),
    winslash = "/",
    mustWork = FALSE
  ))
  if (length(result_files) == 0L) stop("No .RData files found in: ", runs_dir, call. = FALSE)

  true_x <- as.vector(as.matrix(utils::read.table(true_file)))
  rows <- lapply(result_files, function(result_file) {
    tryCatch(
      summarize_mcmc_file(result_file, true_x),
      error = function(err) {
        warning(sprintf("Skipping %s: %s", basename(result_file), conditionMessage(err)), call. = FALSE)
        NULL
      }
    )
  })
  comparison_df <- do.call(rbind, rows[!vapply(rows, is.null, logical(1L))])
  if (is.null(comparison_df) || nrow(comparison_df) == 0L) {
    stop("No valid MCMC runs could be summarized.", call. = FALSE)
  }

  comparison_df <- comparison_df[order(
    comparison_df$rmse_x,
    comparison_df$rmse_beta,
    comparison_df$delta_1,
    comparison_df$delta_2,
    comparison_df$file
  ), ]
  comparison_df$rank_rmse_x <- seq_len(nrow(comparison_df))
  beta_rank <- order(
    comparison_df$rmse_beta,
    comparison_df$rmse_x,
    comparison_df$delta_1,
    comparison_df$delta_2,
    comparison_df$file
  )
  comparison_df$rank_rmse_beta <- NA_integer_
  comparison_df$rank_rmse_beta[beta_rank] <- seq_len(nrow(comparison_df))

  csv_file <- file.path(output_dir, "mcmc_run_comparison_metrics.csv")
  longtable_file <- file.path(output_dir, "mcmc_run_comparison_longtable.tex")
  utils::write.csv(comparison_df, csv_file, row.names = FALSE)
  write_comparison_longtable(
    comparison_df,
    longtable_file,
    sprintf("MCMC runs in %s ranked by RMSE(X).", basename(runs_dir)),
    "tab:mcmc-run-comparison"
  )

  message("Wrote metrics CSV: ", csv_file)
  message("Wrote longtable LaTeX summary: ", longtable_file)
  invisible(list(comparison_df = comparison_df, metrics_csv = csv_file, longtable = longtable_file))
}
