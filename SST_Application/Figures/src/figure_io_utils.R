resolve_figure_chain_file <- function(fig_root,
                                      env_names = character(0),
                                      default_filename = file.path("..", "Data", "results", "sst_application_selected_run.RData")) {
  env_names <- unique(c(env_names, "SST_FIGURE_CHAIN_FILE"))
  env_names <- env_names[nzchar(env_names)]

  for (env_name in env_names) {
    env_value <- Sys.getenv(env_name, unset = "")
    if (nzchar(env_value)) {
      return(normalizePath(env_value, winslash = "/", mustWork = FALSE))
    }
  }

  normalizePath(file.path(fig_root, default_filename), winslash = "/", mustWork = FALSE)
}

figure_run_label <- function(chain_file) {
  label <- basename(chain_file)
  label <- sub("(?i)\\.rdata$", "", label, perl = TRUE)
  label <- gsub("[^A-Za-z0-9._-]+", "_", label)
  if (!nzchar(label)) {
    label <- "run"
  }
  label
}

resolve_figure_output_dir <- function(fig_root, chain_file = NULL, run_label = NULL) {
  output_dir_env <- Sys.getenv("SST_FIGURE_OUTPUT_DIR", unset = "")
  if (nzchar(output_dir_env)) {
    return(normalizePath(output_dir_env, winslash = "/", mustWork = FALSE))
  }

  if (!is.null(run_label) && nzchar(run_label)) {
    return(normalizePath(
      file.path(fig_root, "out", run_label),
      winslash = "/",
      mustWork = FALSE
    ))
  }

  if (!is.null(chain_file) && nzchar(chain_file)) {
    return(normalizePath(
      file.path(fig_root, "out", figure_run_label(chain_file)),
      winslash = "/",
      mustWork = FALSE
    ))
  }

  normalizePath(file.path(fig_root, "out"), winslash = "/", mustWork = FALSE)
}

crop_pdf_in_place <- function(pdf_path) {
  cropped_pdf <- tempfile(
    pattern = paste0(tools::file_path_sans_ext(basename(pdf_path)), "_crop_"),
    tmpdir = dirname(pdf_path),
    fileext = ".pdf"
  )
  crop_status <- suppressWarnings(system2("pdfcrop", args = c(pdf_path, cropped_pdf)))
  if (identical(crop_status, 0L) && file.exists(cropped_pdf)) {
    ok_rename <- file.rename(cropped_pdf, pdf_path)
    if (!isTRUE(ok_rename)) {
      stop("pdfcrop succeeded but failed to replace the original PDF.", call. = FALSE)
    }
  } else if (file.exists(cropped_pdf)) {
    unlink(cropped_pdf)
  }
}
