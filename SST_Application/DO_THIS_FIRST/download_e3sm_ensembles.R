#!/usr/bin/env Rscript

# Download CMIP6 E3SM SST files for the E3SM ensemble members used here.
#
# Examples from Rcode/SST_Application:
#   Rscript DO_THIS_FIRST/download_e3sm_ensembles.R
#   Rscript DO_THIS_FIRST/download_e3sm_ensembles.R --members=r1i1p1f1,r2i1p1f1
#   Rscript DO_THIS_FIRST/download_e3sm_ensembles.R --time_range=201001-201412 --version=auto

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
sst_app_dir <- if (basename(script_dir) == "DO_THIS_FIRST") {
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
} else {
  script_dir
}
setwd(sst_app_dir)

default_members <- sprintf("r%di1p1f1", 1:21)

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript DO_THIS_FIRST/download_e3sm_ensembles.R [member ...] [--key=value ...]\n\n",
    "Members:\n",
    "  Defaults to all configured members: r1i1p1f1 through r21i1p1f1.\n",
    "  Pass positional args or --members=id1,id2,... only to download a subset.\n\n",
    "Options (defaults shown):\n",
    "  --out_dir=Data/raw\n",
    "  --base_url=https://esgf-node.ornl.gov/thredds/fileServer/user_pub_work/CMIP6\n",
    "  --activity=CMIP\n",
    "  --institution=E3SM-Project\n",
    "  --model=E3SM-2-0\n",
    "  --experiment=historical\n",
    "  --table=Omon\n",
    "  --variable=tos\n",
    "  --grid=gr\n",
    "  --version=auto\n",
    "  --time_range=201001-201412\n",
    "  --skip_existing=TRUE\n",
    "  --timeout=600\n",
    "  --help\n",
    sep = ""
  )
}

parse_key_value <- function(arg) {
  kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
  if (length(kv) < 2L) stop("Expected --key=value format, got: ", arg)
  key <- kv[1]
  val <- paste(kv[-1], collapse = "=")
  list(key = key, val = val)
}

split_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

as_bool <- function(x) {
  y <- toupper(trimws(x))
  if (y %in% c("TRUE", "T", "1", "YES", "Y")) return(TRUE)
  if (y %in% c("FALSE", "F", "0", "NO", "N")) return(FALSE)
  stop("Invalid boolean value: ", x, " (use TRUE/FALSE)")
}

args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("-h", "--help"))) {
  usage()
  quit(status = 0)
}

opts <- list(
  out_dir = "Data/raw",
  base_url = "https://esgf-node.ornl.gov/thredds/fileServer/user_pub_work/CMIP6",
  activity = "CMIP",
  institution = "E3SM-Project",
  model = "E3SM-2-0",
  experiment = "historical",
  table = "Omon",
  variable = "tos",
  grid = "gr",
  version = "auto",
  time_range = "201001-201412",
  skip_existing = "TRUE",
  timeout = "600",
  members = ""
)

positional_members <- character(0)
for (arg in args) {
  if (startsWith(arg, "--")) {
    kv <- parse_key_value(arg)
    if (!kv$key %in% names(opts)) stop("Unknown option: --", kv$key)
    opts[[kv$key]] <- kv$val
  } else {
    positional_members <- c(positional_members, arg)
  }
}

members <- unique(c(positional_members, split_csv(opts$members)))
if (length(members) == 0L) {
  members <- default_members
}

skip_existing <- as_bool(opts$skip_existing)
timeout_sec <- as.integer(opts$timeout)
if (is.na(timeout_sec) || timeout_sec <= 0L) stop("--timeout must be a positive integer")

old_timeout <- getOption("timeout")
on.exit(options(timeout = old_timeout), add = TRUE)
options(timeout = max(old_timeout, timeout_sec))

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

build_filename <- function(member) {
  sprintf(
    "%s_%s_%s_%s_%s_%s_%s.nc",
    opts$variable, opts$table, opts$model, opts$experiment,
    member, opts$grid, opts$time_range
  )
}

build_url <- function(member, version) {
  paste(
    opts$base_url,
    opts$activity,
    opts$institution,
    opts$model,
    opts$experiment,
    member,
    opts$table,
    opts$variable,
    opts$grid,
    version,
    build_filename(member),
    sep = "/"
  )
}

catalog_base_url <- sub(
  "/thredds/fileServer/",
  "/thredds/catalog/",
  opts$base_url,
  fixed = TRUE
)
if (identical(catalog_base_url, opts$base_url)) {
  stop("base_url must include '/thredds/fileServer/' so catalog lookup can be derived.")
}

build_catalog_url <- function(member) {
  paste(
    catalog_base_url,
    opts$activity,
    opts$institution,
    opts$model,
    opts$experiment,
    member,
    opts$table,
    opts$variable,
    opts$grid,
    "catalog.xml",
    sep = "/"
  )
}

sort_versions_desc <- function(versions) {
  versions <- unique(versions)
  version_num <- suppressWarnings(as.integer(sub("^v", "", versions)))
  if (all(is.na(version_num))) {
    sort(versions, decreasing = TRUE)
  } else {
    versions[order(version_num, decreasing = TRUE, na.last = TRUE)]
  }
}

get_available_versions <- function(member) {
  catalog_url <- build_catalog_url(member)
  lines <- tryCatch(
    readLines(catalog_url, warn = FALSE),
    error = function(e) character(0)
  )
  if (length(lines) == 0L) return(character(0))
  versions <- regmatches(lines, gregexpr("v[0-9]{8}", lines, perl = TRUE))
  versions <- unlist(versions, use.names = FALSE)
  sort_versions_desc(versions)
}

download_once <- function(url, dest) {
  status <- tryCatch(
    utils::download.file(url = url, destfile = dest, mode = "wb", quiet = FALSE),
    error = function(e) e
  )
  !(inherits(status, "error") || (is.numeric(status) && status != 0))
}

results <- data.frame(
  member = members,
  status = rep("PENDING", length(members)),
  version = rep("", length(members)),
  dest = rep("", length(members)),
  url = rep("", length(members)),
  stringsAsFactors = FALSE
)

for (i in seq_along(members)) {
  member <- members[i]
  filename <- build_filename(member)
  dest <- file.path(opts$out_dir, filename)

  results$dest[i] <- dest

  if (skip_existing && file.exists(dest)) {
    message("[", i, "/", length(members), "] Skip existing: ", basename(dest))
    results$status[i] <- "SKIPPED"
    next
  }

  use_auto_version <- tolower(opts$version) == "auto"
  versions_to_try <- if (use_auto_version) get_available_versions(member) else opts$version

  if (length(versions_to_try) == 0L) {
    message("[", i, "/", length(members), "] FAILED ", member, ": no version folders found")
    results$status[i] <- "FAILED"
    next
  }

  message(
    "[", i, "/", length(members), "] Downloading ", member,
    " (version ", versions_to_try[1], if (length(versions_to_try) > 1L) ", fallback enabled" else "", ")"
  )

  got_file <- FALSE
  for (version in versions_to_try) {
    url <- build_url(member, version)
    ok <- download_once(url, dest)
    if (ok) {
      results$status[i] <- "OK"
      results$version[i] <- version
      results$url[i] <- url
      message("  OK: ", dest)
      got_file <- TRUE
      break
    }
    if (file.exists(dest)) file.remove(dest)
    if (!use_auto_version) {
      results$url[i] <- url
      break
    }
  }

  if (!got_file) {
    if (!nzchar(results$url[i])) {
      failed_version <- versions_to_try[length(versions_to_try)]
      results$url[i] <- build_url(member, failed_version)
    }
    results$version[i] <- paste(versions_to_try, collapse = ",")
    results$status[i] <- "FAILED"
    message("  FAILED: ", results$url[i])
  }
}

cat("\nSummary:\n", sep = "")
print(table(results$status), quote = FALSE)

failed <- results$status == "FAILED"
if (any(failed)) {
  cat("\nFailed members:\n", sep = "")
  print(results[failed, c("member", "version", "url")], row.names = FALSE)
  quit(status = 1)
}

cat("\nAll requested downloads completed.\n")
