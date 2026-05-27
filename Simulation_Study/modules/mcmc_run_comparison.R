# MCMC run comparison helpers.
#
# Source this module and call compare_mcmc_runs(), or use the lightweight
# wrapper at ../compare_mcmc_runs.R.

resolve_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) return(NULL)
  normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE)
}

resolve_project_root <- function() {
  script_path <- resolve_script_path()
  start_dirs <- c(
    if (!is.null(script_path)) dirname(script_path) else NA_character_,
    normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
  )
  start_dirs <- unique(start_dirs[!is.na(start_dirs)])

  parent_dirs <- function(path) {
    dirs <- character(0)
    current <- normalizePath(path, winslash = "/", mustWork = FALSE)
    repeat {
      dirs <- c(dirs, current)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
    dirs
  }

  candidates <- unique(unlist(lapply(start_dirs, parent_dirs), use.names = FALSE))
  candidates <- unique(candidates[!is.na(candidates)])

  for (cand in candidates) {
    if (file.exists(file.path(cand, "run_mcmc.R")) &&
        file.exists(file.path(cand, "modules", "spectral2_ens.R")) &&
        dir.exists(file.path(cand, "Data", "inputs"))) {
      return(cand)
    }
  }

  stop("Could not locate Simulation_Study root containing run_mcmc.R, modules/spectral2_ens.R, and Data/inputs/.", call. = FALSE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

configure_project_lib <- function(project_root) {
  project_lib <- file.path(project_root, ".Rlib")
  if (dir.exists(project_lib)) {
    .libPaths(c(project_lib, .libPaths()))
  }
}

load_comparison_packages <- function() {
  suppressPackageStartupMessages({
    library(ggplot2)
  })
}

is_absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[\\\\/])", path)
}

resolve_path <- function(path, base_dir) {
  path <- path.expand(path)
  if (is_absolute_path(path)) {
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
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

safe_rmse <- function(truth, estimate) {
  truth <- as.numeric(truth)
  estimate <- as.numeric(estimate)
  ok <- is.finite(truth) & is.finite(estimate)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((estimate[ok] - truth[ok])^2))
}

read_true_field <- function(path) {
  true_mat <- as.matrix(utils::read.table(path))
  as.vector(true_mat)
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

extract_timestamp <- function(file_name) {
  m <- regmatches(file_name, regexpr("[0-9]{8}_[0-9]{6}_[0-9]{6}", file_name))
  if (length(m) == 0L || !nzchar(m)) return(sub("\\.RData$", "", file_name))
  m
}

short_run_label <- function(file_name, setup, params, settings) {
  k0 <- as.integer(setup$K0 %||% setting_int(settings, "k0"))
  seed <- as.integer(setup$seed %||% setting_int(settings, "seed"))
  pk <- as.numeric(unique(params$p_k %||% setup$p_k_const %||% setting_num(settings, "p_k"))[1L])
  d1 <- as.numeric(params$delta1 %||% setting_num(settings, "delta1"))
  d2 <- as.numeric(params$delta2 %||% setting_num(settings, "delta2"))

  pieces <- c(
    if (is.finite(k0)) sprintf("K0=%d", k0) else NULL,
    if (is.finite(d1)) sprintf("d1=%s", format_decimal(d1, 2L)) else NULL,
    if (is.finite(d2)) sprintf("d2=%s", format_decimal(d2, 2L)) else NULL,
    if (is.finite(pk)) sprintf("p=%s", format_decimal(pk, 2L)) else NULL,
    if (is.finite(seed)) sprintf("seed=%d", seed) else NULL
  )

  if (length(pieces) == 0L) {
    return(extract_timestamp(file_name))
  }
  paste(pieces, collapse = ", ")
}

extract_run_summary <- function(result_file, true_full) {
  settings <- read_settings_file(result_file)

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

  Yc_bar <- as.numeric(static$Yc_bar %||% numeric(0))
  if (length(Yc_bar) == 0L) stop("missing static$Yc_bar")

  delta1 <- as.numeric(params$delta1 %||% setting_num(settings, "delta1"))
  delta2 <- as.numeric(params$delta2 %||% setting_num(settings, "delta2"))
  k0 <- as.integer(setup$K0 %||% setting_int(settings, "k0"))
  seed <- as.integer(setup$seed %||% setting_int(settings, "seed"))
  pk <- as.numeric(unique(params$p_k %||% setup$p_k_const %||% setting_num(settings, "p_k"))[1L])
  m_ens <- as.integer(setup$ensemble_count %||% setting_int(settings, "m_ens"))
  burn_in <- as.integer(run_config$burn_in %||% 0L)
  update_sigma2_after <- as.integer(run_config$update_sigma2_after %||% burn_in)
  summary_after <- max(c(0L, burn_in, update_sigma2_after), na.rm = TRUE)

  n_sites <- length(true_full)
  if (length(Yc_bar) != n_sites) {
    stop(sprintf("Yc_bar length %d does not match trueXX length %d", length(Yc_bar), n_sites))
  }

  iter_chain <- as.integer(store$iter_chain %||% integer(0))
  keep_idx <- as.integer(output$keep_idx %||% seq_along(iter_chain))
  keep_idx <- unique(keep_idx[is.finite(keep_idx) & keep_idx >= 1L])
  if (length(iter_chain) > 0L) {
    keep_idx <- keep_idx[keep_idx <= length(iter_chain)]
  }
  if (length(keep_idx) == 0L) stop("no kept MCMC draws available")

  comparison_keep_idx <- keep_idx
  if (length(iter_chain) > 0L) {
    comparison_keep_idx <- keep_idx[iter_chain[keep_idx] > summary_after]
    if (length(comparison_keep_idx) == 0L) {
      comparison_keep_idx <- tail(keep_idx, 1L)
    }
  }

  chain_post_mean <- function(chain_matrix, idx, fallback) {
    chain_matrix <- as.matrix(chain_matrix %||% matrix(numeric(0), nrow = 0L, ncol = 0L))
    idx <- idx[idx <= ncol(chain_matrix)]
    if (length(idx) == 0L) {
      return(as.numeric(fallback %||% numeric(0)))
    }
    rowMeans(chain_matrix[, idx, drop = FALSE], na.rm = TRUE)
  }

  kept_mean <- function(x) {
    x <- as.numeric(x %||% numeric(0))
    idx <- comparison_keep_idx[comparison_keep_idx <= length(x)]
    if (length(idx) == 0L) return(NA_real_)
    mean(x[idx], na.rm = TRUE)
  }

  X_post_mean <- chain_post_mean(store$X_chain, comparison_keep_idx, output$X_post_mean)
  beta_post_mean <- chain_post_mean(store$beta_chain, comparison_keep_idx, output$beta_post_mean)

  if (length(X_post_mean) == 0L) stop("missing output/store X posterior draws")
  if (length(beta_post_mean) == 0L) stop("missing output/store beta posterior draws")
  if (length(X_post_mean) != n_sites) {
    stop(sprintf("X_post_mean length %d does not match trueXX length %d", length(X_post_mean), n_sites))
  }
  if (length(beta_post_mean) != n_sites) {
    stop(sprintf("beta_post_mean length %d does not match trueXX length %d", length(beta_post_mean), n_sites))
  }

  true_beta <- Yc_bar - true_full
  file_name <- basename(result_file)
  run_id <- extract_timestamp(file_name)

  data.frame(
    file = file_name,
    run = run_id,
    run_label = short_run_label(file_name, setup, params, settings),
    timestamp = run_id,
    K0 = k0,
    M_ens = m_ens,
    p_k = pk,
    seed = seed,
    delta_1 = delta1,
    delta_2 = delta2,
    rmse_x = safe_rmse(true_full, X_post_mean),
    rmse_beta = safe_rmse(true_beta, beta_post_mean),
    post_mean_sigma2yc = kept_mean(store$sig2c_chain),
    post_mean_sigma2yo = kept_mean(store$sig2o_chain),
    post_mean_d_beta = kept_mean(store$d_beta_chain),
    n_kept = length(comparison_keep_idx),
    sig2yc_init = as.numeric(run_config$sig2c_init %||% setting_num(settings, "SIG2C")),
    sig2yo_init = as.numeric(run_config$sig2o_init %||% setting_num(settings, "SIG2O")),
    stringsAsFactors = FALSE
  )
}

compare_mcmc_runs <- function(runs_dir = NULL,
                              output_dir = NULL,
                              true_file = NULL,
                              invocation_dir = normalizePath(getwd(), winslash = "/", mustWork = FALSE),
                              project_root = NULL) {
  invocation_dir <- normalizePath(invocation_dir, winslash = "/", mustWork = FALSE)
  project_root <- project_root %||% resolve_project_root()

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)

  configure_project_lib(project_root)
  load_comparison_packages()

  runs_dir <- if (is.null(runs_dir) || length(runs_dir) == 0L || !nzchar(runs_dir[1L])) {
    file.path(project_root, "Data", "results")
  } else {
    resolve_path(runs_dir[1L], invocation_dir)
  }
  runs_dir <- normalizePath(runs_dir, winslash = "/", mustWork = FALSE)

  output_dir <- if (is.null(output_dir) || length(output_dir) == 0L || !nzchar(output_dir[1L])) {
    file.path(runs_dir, "comparison")
  } else {
    resolve_path(output_dir[1L], invocation_dir)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  true_file <- if (is.null(true_file) || length(true_file) == 0L || !nzchar(true_file[1L])) {
    file.path(project_root, "Data", "inputs", "trueXX.txt")
  } else {
    resolve_path(true_file[1L], invocation_dir)
  }
  true_file <- normalizePath(true_file, winslash = "/", mustWork = FALSE)

  if (!dir.exists(runs_dir)) {
    stop("MCMC runs directory does not exist: ", runs_dir, call. = FALSE)
  }
  if (!file.exists(true_file)) {
    stop("trueXX file does not exist: ", true_file, call. = FALSE)
  }

  result_files <- list.files(
    runs_dir,
    pattern = "\\.RData$",
    recursive = FALSE,
    full.names = TRUE
  )
  result_files <- sort(normalizePath(result_files, winslash = "/", mustWork = FALSE))

  if (length(result_files) == 0L) {
    stop("No .RData files found in: ", runs_dir, call. = FALSE)
  }

  true_full <- read_true_field(true_file)
  
  message("Comparing ", length(result_files), " MCMC result files in: ", runs_dir)
  rows <- lapply(result_files, function(result_file) {
    tryCatch(
      extract_run_summary(result_file, true_full = true_full),
      error = function(err) {
        warning(
          sprintf("Skipping %s: %s", basename(result_file), conditionMessage(err)),
          call. = FALSE
        )
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
  x_table_file <- file.path(output_dir, "top10_rmse_x_table.tex")
  x_table_with_run_file <- file.path(output_dir, "top10_rmse_x_table_with_run.tex")
  x_table_without_run_file <- file.path(output_dir, "top10_rmse_x_table_without_run.tex")
  beta_table_file <- file.path(output_dir, "top10_rmse_beta_table.tex")
  beta_table_with_run_file <- file.path(output_dir, "top10_rmse_beta_table_with_run.tex")
  beta_table_without_run_file <- file.path(output_dir, "top10_rmse_beta_table_without_run.tex")
  tables_file <- file.path(output_dir, "top10_rmse_tables.tex")
  tables_with_run_file <- file.path(output_dir, "top10_rmse_tables_with_run.tex")
  tables_without_run_file <- file.path(output_dir, "top10_rmse_tables_without_run.tex")
  x_figure_pdf <- file.path(output_dir, "top10_rmse_x_plot.pdf")
  x_figure_png <- file.path(output_dir, "top10_rmse_x_plot.png")
  beta_figure_pdf <- file.path(output_dir, "top10_rmse_beta_plot.pdf")
  beta_figure_png <- file.path(output_dir, "top10_rmse_beta_plot.png")
  figure_pdf <- file.path(output_dir, "top10_rmse_4panel_plot.pdf")
  figure_png <- file.path(output_dir, "top10_rmse_4panel_plot.png")
  longtable_file <- file.path(output_dir, "mcmc_run_comparison_longtable.tex")
  tradeoff_pdf <- file.path(output_dir, "rmse_tradeoff_plot.pdf")
  tradeoff_png <- file.path(output_dir, "rmse_tradeoff_plot.png")
  
  utils::write.csv(comparison_df, csv_file, row.names = FALSE)
  
  latex_bold <- function(x) {
    paste0("\\textbf{", x, "}")
  }
  
  make_top10_table_lines <- function(top_runs, rank_col, common_runs, title_comment, include_run = TRUE) {
    rank_values <- top_runs[[rank_col]]
  
    if (any(top_runs$run %in% common_runs)) {
      common_note <- "% Common runs between top-10 RMSE(X) and top-10 RMSE(beta) are bolded."
    } else {
      common_note <- "% No common runs between top-10 RMSE(X) and top-10 RMSE(beta)."
    }
  
    if (include_run) {
      header <- c(
        "Rank",
        "Run",
        "RMSE$(X)$",
        "RMSE$(\\beta)$",
        "$\\delta_1$",
        "$\\delta_2$",
        "$\\overline{\\sigma^2_{yc}}$",
        "$\\overline{\\sigma^2_{yo}}$",
        "$\\overline{d_\\beta}$"
      )
      table_lines <- c(
        title_comment,
        common_note,
        "\\begin{tabular}{@{}rlrrrrrrr@{}}",
        "\\toprule",
        paste0(paste(header, collapse = " & "), " \\\\"),
        "\\midrule"
      )
    } else {
      header <- c(
        "Rank",
        "RMSE$(X)$",
        "RMSE$(\\beta)$",
        "$\\delta_1$",
        "$\\delta_2$",
        "$\\overline{\\sigma^2_{yc}}$",
        "$\\overline{\\sigma^2_{yo}}$",
        "$\\overline{d_\\beta}$"
      )
      table_lines <- c(
        title_comment,
        common_note,
        "\\begin{tabular}{@{}rrrrrrrr@{}}",
        "\\toprule",
        paste0(paste(header, collapse = " & "), " \\\\"),
        "\\midrule"
      )
    }
  
    for (i in seq_len(nrow(top_runs))) {
      row <- top_runs[i, ]
      cells <- c(
        as.character(rank_values[i]),
        if (include_run) latex_escape(row$run) else NULL,
        format_decimal(row$rmse_x, 3L),
        format_decimal(row$rmse_beta, 3L),
        format_decimal(row$delta_1, 3L),
        format_decimal(row$delta_2, 3L),
        format_decimal(row$post_mean_sigma2yc, 3L),
        format_decimal(row$post_mean_sigma2yo, 3L),
        format_decimal(row$post_mean_d_beta, 2L)
      )
      if (row$run %in% common_runs) {
        cells <- latex_bold(cells)
      }
      table_lines <- c(table_lines, paste0(paste(cells, collapse = " & "), " \\\\"))
    }
  
    c(
      table_lines,
      "\\bottomrule",
      "\\end{tabular}"
    )
  }
  
  write_top10_tables <- function(comparison_df) {
    top_x <- comparison_df[order(comparison_df$rank_rmse_x), , drop = FALSE]
    top_x <- head(top_x, 10L)
    top_beta <- comparison_df[order(comparison_df$rank_rmse_beta), , drop = FALSE]
    top_beta <- head(top_beta, 10L)
    common_runs <- intersect(top_x$run, top_beta$run)
  
    x_table_lines_with_run <- make_top10_table_lines(
      top_runs = top_x,
      rank_col = "rank_rmse_x",
      common_runs = common_runs,
      title_comment = "% Top 10 runs by RMSE(X), with Run column",
      include_run = TRUE
    )
    beta_table_lines_with_run <- make_top10_table_lines(
      top_runs = top_beta,
      rank_col = "rank_rmse_beta",
      common_runs = common_runs,
      title_comment = "% Top 10 runs by RMSE(beta), with Run column",
      include_run = TRUE
    )
    x_table_lines_without_run <- make_top10_table_lines(
      top_runs = top_x,
      rank_col = "rank_rmse_x",
      common_runs = common_runs,
      title_comment = "% Top 10 runs by RMSE(X), without Run column",
      include_run = FALSE
    )
    beta_table_lines_without_run <- make_top10_table_lines(
      top_runs = top_beta,
      rank_col = "rank_rmse_beta",
      common_runs = common_runs,
      title_comment = "% Top 10 runs by RMSE(beta), without Run column",
      include_run = FALSE
    )
  
    writeLines(x_table_lines_with_run, x_table_file)
    writeLines(x_table_lines_with_run, x_table_with_run_file)
    writeLines(x_table_lines_without_run, x_table_without_run_file)
    writeLines(beta_table_lines_with_run, beta_table_file)
    writeLines(beta_table_lines_with_run, beta_table_with_run_file)
    writeLines(beta_table_lines_without_run, beta_table_without_run_file)
    writeLines(c(x_table_lines_with_run, "", beta_table_lines_with_run), tables_file)
    writeLines(c(x_table_lines_with_run, "", beta_table_lines_with_run), tables_with_run_file)
    writeLines(c(x_table_lines_without_run, "", beta_table_lines_without_run), tables_without_run_file)
  
    invisible(common_runs)
  }
  
  make_top10_plots <- function(comparison_df, output_x_pdf, output_x_png, output_beta_pdf, output_beta_png, output_four_pdf, output_four_png) {
    top_x <- comparison_df[order(comparison_df$rank_rmse_x), , drop = FALSE]
    top_x <- head(top_x, 10L)
    top_x$x_rank_label <- sprintf("%02d", top_x$rank_rmse_x)
    top_x$plot_label <- top_x$x_rank_label
    top_x_axis_labels <- parse(text = sprintf(
      "paste('%s', '  ', bar(sigma[yo]^2), ' = %s')",
      top_x$x_rank_label,
      format_decimal(top_x$post_mean_sigma2yo, 2L)
    ))
    names(top_x_axis_labels) <- top_x$plot_label
    top_x$plot_label <- factor(top_x$plot_label, levels = rev(top_x$plot_label))
    x_min_rmse_x <- min(top_x$rmse_x, na.rm = TRUE) * 0.995
  
    rank_plot <- ggplot(top_x, aes(x = rmse_x, y = plot_label)) +
      geom_segment(
        aes(xend = rmse_x, yend = plot_label),
        x = x_min_rmse_x,
        linewidth = 0.45,
        color = "grey65"
      ) +
      geom_point(aes(color = rmse_beta, size = post_mean_d_beta), alpha = 0.92) +
      scale_color_gradient(
        low = "#2c7fb8",
        high = "#d95f02",
        name = expression("RMSE(" * beta * ")")
      ) +
      scale_size_continuous(
        range = c(2.8, 7.2),
        name = expression(bar(d)[beta])
      ) +
      scale_x_continuous(expand = expansion(mult = c(0.01, 0.08))) +
      scale_y_discrete(labels = top_x_axis_labels) +
      labs(
        title = "(A) Top 10 runs by X RMSE",
        x = "RMSE(X)",
        y = NULL
      ) +
      theme_bw(base_size = 9.5, base_family = "serif") +
      theme(
        legend.position = "right",
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7.4, color = "black"),
        axis.text.x = element_text(color = "black"),
        plot.title = element_text(face = "bold", hjust = 0),
        plot.margin = margin(5.5, 5.5, 5.5, 18)
      )
  
    tradeoff_plot <- ggplot(top_x, aes(x = rmse_x, y = rmse_beta)) +
      geom_path(
        aes(group = 1L),
        color = "grey70",
        linewidth = 0.35,
        linetype = "dotted"
      ) +
      geom_point(
        color = "#1b9e77",
        fill = "#1b9e77",
        size = 6.2,
        alpha = 0.88
      ) +
      geom_text(
        aes(label = x_rank_label),
        color = "white",
        fontface = "bold",
        size = 3.0
      ) +
      labs(
        title = "(B) Beta RMSE among X-best runs",
        x = "RMSE(X)",
        y = expression("RMSE(" * beta * ")")
      ) +
      theme_bw(base_size = 9.5, base_family = "serif") +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", hjust = 0)
      )
  
    figure_grob <- gridExtra::arrangeGrob(
      rank_plot,
      tradeoff_plot,
      nrow = 1,
      widths = c(1.35, 1)
    )
  
    top_beta <- comparison_df[order(comparison_df$rank_rmse_beta), , drop = FALSE]
    top_beta <- head(top_beta, 10L)
    top_beta$beta_rank_label <- sprintf("%02d", top_beta$rank_rmse_beta)
    top_beta$plot_label <- top_beta$beta_rank_label
    top_beta_axis_labels <- parse(text = sprintf(
      "paste('%s', '  ', bar(sigma[yo]^2), ' = %s')",
      top_beta$beta_rank_label,
      format_decimal(top_beta$post_mean_sigma2yo, 2L)
    ))
    names(top_beta_axis_labels) <- top_beta$plot_label
    top_beta$plot_label <- factor(top_beta$plot_label, levels = rev(top_beta$plot_label))
    x_min_rmse_beta <- min(top_beta$rmse_beta, na.rm = TRUE) * 0.995
  
    beta_rank_plot <- ggplot(top_beta, aes(x = rmse_beta, y = plot_label)) +
      geom_segment(
        aes(xend = rmse_beta, yend = plot_label),
        x = x_min_rmse_beta,
        linewidth = 0.45,
        color = "grey65"
      ) +
      geom_point(aes(color = rmse_x, size = post_mean_d_beta), alpha = 0.92) +
      scale_color_gradient(
        low = "#2c7fb8",
        high = "#d95f02",
        name = "RMSE(X)"
      ) +
      scale_size_continuous(
        range = c(2.8, 7.2),
        name = expression(bar(d)[beta])
      ) +
      scale_x_continuous(expand = expansion(mult = c(0.01, 0.08))) +
      scale_y_discrete(labels = top_beta_axis_labels) +
      labs(
        title = "(A) Top 10 runs by beta RMSE",
        x = expression("RMSE(" * beta * ")"),
        y = NULL
      ) +
      theme_bw(base_size = 9.5, base_family = "serif") +
      theme(
        legend.position = "right",
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7.4, color = "black"),
        axis.text.x = element_text(color = "black"),
        plot.title = element_text(face = "bold", hjust = 0),
        plot.margin = margin(5.5, 5.5, 5.5, 18)
      )
  
    beta_tradeoff_plot <- ggplot(top_beta, aes(x = rmse_beta, y = rmse_x)) +
      geom_path(
        aes(group = 1L),
        color = "grey70",
        linewidth = 0.35,
        linetype = "dotted"
      ) +
      geom_point(
        color = "#1b9e77",
        fill = "#1b9e77",
        size = 6.2,
        alpha = 0.88
      ) +
      geom_text(
        aes(label = beta_rank_label),
        color = "white",
        fontface = "bold",
        size = 3.0
      ) +
      labs(
        title = "(B) X RMSE among beta-best runs",
        x = expression("RMSE(" * beta * ")"),
        y = "RMSE(X)"
      ) +
      theme_bw(base_size = 9.5, base_family = "serif") +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", hjust = 0)
      )
  
    beta_figure_grob <- gridExtra::arrangeGrob(
      beta_rank_plot,
      beta_tradeoff_plot,
      nrow = 1,
      widths = c(1.35, 1)
    )
  
    compact_four_panel_theme <- theme(
      axis.text = element_text(color = "black", size = 7.0),
      axis.text.y = element_text(color = "black", size = 6.5),
      axis.title = element_text(size = 8.2),
      plot.title = element_text(face = "bold", hjust = 0, size = 9.2),
      plot.margin = margin(2.5, 2.5, 2.5, 7),
      panel.grid.minor = element_blank()
    )
  
    compact_lollipop_guides <- guides(
      color = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(0.20, "in"),
        barheight = grid::unit(0.90, "in")
      ),
      size = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        keyheight = grid::unit(0.14, "in"),
        keywidth = grid::unit(0.22, "in"),
        override.aes = list(alpha = 0.82)
      )
    )
  
    compact_lollipop_legend_theme <- theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = grid::unit(0.08, "in"),
      legend.key.size = grid::unit(0.16, "in"),
      legend.title = element_text(size = 6.6),
      legend.text = element_text(size = 6.1)
    )
  
    compact_no_legend_theme <- theme(legend.position = "none")
  
    figure_4panel_grob <- gridExtra::arrangeGrob(
      rank_plot +
        labs(title = "(A) Ranked by RMSE(X)") +
        compact_lollipop_guides +
        compact_four_panel_theme +
        compact_lollipop_legend_theme,
      tradeoff_plot +
        labs(title = "(B) X-best tradeoff") +
        compact_four_panel_theme +
        compact_no_legend_theme,
      beta_rank_plot +
        labs(title = "(C) Ranked by RMSE(beta)") +
        compact_lollipop_guides +
        compact_four_panel_theme +
        compact_lollipop_legend_theme,
      beta_tradeoff_plot +
        labs(title = "(D) Beta-best tradeoff") +
        compact_four_panel_theme +
        compact_no_legend_theme,
      nrow = 2,
      ncol = 2,
      widths = c(1.55, 1.43),
      heights = c(1, 1),
      padding = grid::unit(0, "pt")
    )
  
    pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggplot2::ggsave(
      filename = output_x_pdf,
      plot = figure_grob,
      device = pdf_device,
      width = 11.0,
      height = 4.6,
      units = "in",
      bg = "white"
    )
    ggplot2::ggsave(
      filename = output_x_png,
      plot = figure_grob,
      width = 11.0,
      height = 4.6,
      units = "in",
      dpi = 400,
      bg = "white"
    )
    ggplot2::ggsave(
      filename = output_beta_pdf,
      plot = beta_figure_grob,
      device = pdf_device,
      width = 11.0,
      height = 4.6,
      units = "in",
      bg = "white"
    )
    ggplot2::ggsave(
      filename = output_beta_png,
      plot = beta_figure_grob,
      width = 11.0,
      height = 4.6,
      units = "in",
      dpi = 400,
      bg = "white"
    )
    ggplot2::ggsave(
      filename = output_four_pdf,
      plot = figure_4panel_grob,
      device = pdf_device,
      width = 10.8,
      height = 6.9,
      units = "in",
      bg = "white"
    )
    ggplot2::ggsave(
      filename = output_four_png,
      plot = figure_4panel_grob,
      width = 10.8,
      height = 6.9,
      units = "in",
      dpi = 400,
      bg = "white"
    )
  }
  
  write_longtable <- function(comparison_df, output_file, caption, label) {
    header_cells <- c(
      "Rank",
      "RMSE$(X)$",
      "RMSE$(\\beta)$",
      "Rank$(\\beta)$",
      "$\\overline{\\sigma^2_{yo}}$",
      "$\\overline{\\sigma^2_{yc}}$",
      "$p_k$",
      "$\\overline{d_{\\beta}}$",
      "$\\delta_1$",
      "$\\delta_2$"
    )
    header_line <- paste0(paste(header_cells, collapse = " & "), " \\\\")
  
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
  
    lines <- c(
      "% Requires \\usepackage{booktabs,longtable}",
      "% Generated by compare_mcmc_runs.R",
      "\\begin{longtable}{@{}rrrrrrrrrr@{}}",
      paste0("\\caption{", latex_escape(caption), "}\\label{", label, "}\\\\"),
      "\\toprule",
      header_line,
      "\\midrule",
      "\\endfirsthead",
      "\\multicolumn{10}{l}{\\tablename\\ \\thetable\\ -- continued from previous page}\\\\",
      "\\toprule",
      header_line,
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
    )
  
    writeLines(lines, output_file)
  }
  
  make_tradeoff_plot <- function(comparison_df, output_pdf, output_png) {
    pk_levels <- sort(unique(comparison_df$p_k[is.finite(comparison_df$p_k)]))
    pk_labels <- format_decimal(pk_levels, 2L)
    comparison_df$pk_label <- factor(
      format_decimal(comparison_df$p_k, 2L),
      levels = pk_labels
    )
    comparison_df$avg_rank <- (comparison_df$rank_rmse_x + comparison_df$rank_rmse_beta) / 2
    comparison_df$delta_label <- sprintf(
      "(%s, %s)",
      format_decimal(comparison_df$delta_1, 2L),
      format_decimal(comparison_df$delta_2, 2L)
    )
  
    if (length(pk_labels) == 0L) {
      comparison_df$pk_label <- factor("all", levels = "all")
      pk_labels <- "all"
    }
  
    facet_highlights <- lapply(split(comparison_df, comparison_df$pk_label), function(df) {
      df <- df[order(df$avg_rank, df$rmse_x, df$rmse_beta), , drop = FALSE]
      df <- head(df, 2L)
      df$facet_rank <- seq_len(nrow(df))
      df
    })
    highlight_df <- do.call(rbind, facet_highlights)
    row.names(highlight_df) <- NULL
    highlight_df$label_x <- highlight_df$rmse_x + c(-0.18, 0.18)[highlight_df$facet_rank]
    highlight_df$label_y <- highlight_df$rmse_beta + c(0.12, -0.12)[highlight_df$facet_rank]
  
    best_run <- comparison_df[comparison_df$rank_rmse_x == 1L, , drop = FALSE]
    x_limits <- range(comparison_df$rmse_x, na.rm = TRUE)
    y_limits <- range(comparison_df$rmse_beta, na.rm = TRUE)
    lim_min <- floor(min(x_limits[1L], y_limits[1L]) * 10) / 10
    lim_max <- ceiling(max(x_limits[2L], y_limits[2L]) * 10) / 10
    panel_mid <- (lim_min + lim_max) / 2
    x_dir <- ifelse(highlight_df$rmse_x <= panel_mid, 1, -1)
    y_dir <- ifelse(highlight_df$rmse_beta <= panel_mid, 1, -1)
    highlight_df$label_x <- highlight_df$rmse_x + x_dir * c(0.30, 0.48)[highlight_df$facet_rank]
    highlight_df$label_y <- highlight_df$rmse_beta + y_dir * c(0.16, 0.30)[highlight_df$facet_rank]
    highlight_df$label_x <- pmin(pmax(highlight_df$label_x, lim_min + 0.12), lim_max - 0.12)
    highlight_df$label_y <- pmin(pmax(highlight_df$label_y, lim_min + 0.12), lim_max - 0.12)
  
    plot_obj <- ggplot(comparison_df, aes(x = rmse_x, y = rmse_beta)) +
      geom_abline(
        intercept = 0,
        slope = 1,
        color = "grey78",
        linewidth = 0.45,
        linetype = "dashed"
      ) +
      geom_point(
        aes(fill = delta_1, size = delta_2),
        shape = 21,
        color = "grey15",
        alpha = 0.88,
        stroke = 0.28
      ) +
      geom_segment(
        data = highlight_df,
        aes(x = rmse_x, y = rmse_beta, xend = label_x, yend = label_y),
        inherit.aes = FALSE,
        color = "grey55",
        linewidth = 0.26
      ) +
      geom_label(
        data = highlight_df,
        aes(x = label_x, y = label_y, label = delta_label),
        inherit.aes = FALSE,
        size = 3.0,
        label.size = 0.14,
        fill = grDevices::adjustcolor("white", alpha.f = 0.92),
        color = "black",
        label.padding = grid::unit(0.10, "lines")
      ) +
      geom_point(
        data = best_run,
        aes(x = rmse_x, y = rmse_beta),
        inherit.aes = FALSE,
        shape = 8,
        size = 4.4,
        stroke = 0.9,
        color = "#c98a00"
      ) +
      scale_fill_gradientn(
        colours = c("#0b4f6c", "#3c8d40", "#f1b722", "#c24e00"),
        breaks = sort(unique(comparison_df$delta_1)),
        name = expression(delta[1])
      ) +
      scale_size_continuous(
        range = c(2.8, 7.0),
        breaks = sort(unique(comparison_df$delta_2)),
        name = expression(delta[2])
      ) +
      scale_x_continuous(
        limits = c(lim_min, lim_max),
        expand = expansion(mult = c(0.02, 0.04))
      ) +
      scale_y_continuous(
        limits = c(lim_min, lim_max),
        expand = expansion(mult = c(0.02, 0.06))
      ) +
      coord_equal() +
      facet_wrap(
        ~ pk_label,
        nrow = 1,
        labeller = ggplot2::as_labeller(
          stats::setNames(
            paste0("p_k = ", pk_labels),
            pk_labels
          )
        )
      ) +
      labs(
        title = "MCMC RMSE Tradeoff by p_k",
        subtitle = paste0(
          "Facets separate p_k to reduce overplotting. Fill is delta_1, size is delta_2, ",
          "the dashed line is RMSE(beta) = RMSE(X), and labels mark the two best average-rank settings in each panel."
        ),
        x = "RMSE(X)",
        y = expression("RMSE(" * beta * ")")
      ) +
      theme_bw(base_size = 10.5, base_family = "serif") +
      theme(
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9.0),
        strip.background = element_rect(fill = "grey97", color = "grey82"),
        strip.text = element_text(face = "bold"),
        plot.margin = margin(8, 16, 8, 8)
      )
  
    pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggplot2::ggsave(
      filename = output_pdf,
      plot = plot_obj,
      device = pdf_device,
      width = 10.2,
      height = 5.4,
      units = "in",
      bg = "white"
    )
    ggplot2::ggsave(
      filename = output_png,
      plot = plot_obj,
      width = 10.2,
      height = 5.4,
      units = "in",
      dpi = 400,
      bg = "white"
    )
  }
  
  common_top10_runs <- write_top10_tables(comparison_df)
  make_top10_plots(
    comparison_df = comparison_df,
    output_x_pdf = x_figure_pdf,
    output_x_png = x_figure_png,
    output_beta_pdf = beta_figure_pdf,
    output_beta_png = beta_figure_png,
    output_four_pdf = figure_pdf,
    output_four_png = figure_png
  )
  
  caption_text <- sprintf("MCMC runs in %s ranked by RMSE(X).", basename(runs_dir))
  write_longtable(
    comparison_df = comparison_df,
    output_file = longtable_file,
    caption = caption_text,
    label = "tab:mcmc-run-comparison"
  )
  make_tradeoff_plot(
    comparison_df = comparison_df,
    output_pdf = tradeoff_pdf,
    output_png = tradeoff_png
  )
  
  message("Wrote metrics CSV: ", csv_file)
  message("Wrote X-RMSE LaTeX table: ", x_table_file)
  message("Wrote X-RMSE LaTeX table with Run column: ", x_table_with_run_file)
  message("Wrote X-RMSE LaTeX table without Run column: ", x_table_without_run_file)
  message("Wrote beta-RMSE LaTeX table: ", beta_table_file)
  message("Wrote beta-RMSE LaTeX table with Run column: ", beta_table_with_run_file)
  message("Wrote beta-RMSE LaTeX table without Run column: ", beta_table_without_run_file)
  message("Wrote combined LaTeX tables: ", tables_file)
  message("Wrote combined LaTeX tables with Run column: ", tables_with_run_file)
  message("Wrote combined LaTeX tables without Run column: ", tables_without_run_file)
  message("Common top-10 runs bolded in tables: ", paste(common_top10_runs, collapse = ", "))
  message("Wrote X-RMSE figure PDF: ", x_figure_pdf)
  message("Wrote X-RMSE figure PNG: ", x_figure_png)
  message("Wrote beta-RMSE figure PDF: ", beta_figure_pdf)
  message("Wrote beta-RMSE figure PNG: ", beta_figure_png)
  message("Wrote 4-panel figure PDF: ", figure_pdf)
  message("Wrote 4-panel figure PNG: ", figure_png)
  message("Wrote longtable LaTeX summary: ", longtable_file)
  message("Wrote RMSE tradeoff figure PDF: ", tradeoff_pdf)
  message("Wrote RMSE tradeoff figure PNG: ", tradeoff_png)
  
  invisible(list(
    comparison_df = comparison_df,
    output_dir = output_dir,
    metrics_csv = csv_file,
    top10_x_table = x_table_file,
    top10_beta_table = beta_table_file,
    longtable = longtable_file,
    tradeoff_pdf = tradeoff_pdf,
    tradeoff_png = tradeoff_png
  ))
}

run_compare_mcmc_runs_cli <- function(args = commandArgs(trailingOnly = TRUE),
                                      invocation_dir = normalizePath(getwd(), winslash = "/", mustWork = FALSE)) {
  project_root <- resolve_project_root()
  compare_mcmc_runs(
    runs_dir = if (length(args) >= 1L) args[1L] else NULL,
    output_dir = if (length(args) >= 2L) args[2L] else NULL,
    true_file = if (length(args) >= 3L) args[3L] else NULL,
    invocation_dir = invocation_dir,
    project_root = project_root
  )
}
