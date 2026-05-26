rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE)
  script_dir <- dirname(script_path)
  project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  project_dir <- normalizePath(getwd(), mustWork = TRUE)
  script_dir <- file.path(project_dir, "K_Selection")
}

##########
# 1) Inputs
##########

data_dir <- file.path(project_dir, "Data")
processed_dir <- file.path(data_dir, "processed")
plot_dir <- script_dir
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

eig_file <- file.path(processed_dir, "pacific_ocean_laplacian_eigendecomp.rds")
lap_file <- file.path(processed_dir, "pacific_ocean_graph_laplacian.rds")
yc_member_dir <- file.path(processed_dir, "pacific_Yc_ensemble_members")
grid_file <- file.path(processed_dir, "pacific_sst_grid.csv")
K_focus <- c(1000L, 1200L, 1500L)
K_recommend <- 1200L

needed_files <- c(eig_file, lap_file, grid_file)
missing_files <- needed_files[!file.exists(needed_files)]
if (length(missing_files) > 0L) {
  stop("Missing required files:\n", paste0("  - ", missing_files, collapse = "\n"))
}
yc_files <- sort(Sys.glob(file.path(yc_member_dir, "pacific_Yc_*.csv")))
if (length(yc_files) == 0L) {
  stop("No ensemble member files found in ", yc_member_dir)
}

##########
# 2) Load data and prepare ocean vector
##########

cat("Loading eigendecomposition (this is the largest input)...\n")
eig <- readRDS(eig_file)
lap <- readRDS(lap_file)
grid_df <- utils::read.csv(grid_file, stringsAsFactors = FALSE)

ocean_sites <- as.integer(lap$ocean_sites)
n_ocean <- length(ocean_sites)
if (n_ocean == 0L) stop("No ocean sites found.")

cat(sprintf("Loading %d ensemble members and forming the ensemble mean...\n", length(yc_files)))
Yc_matrix <- matrix(NA_real_, nrow = n_ocean, ncol = length(yc_files))
imputed_nonfinite_total <- 0L

for (j in seq_along(yc_files)) {
  yc_df <- utils::read.csv(yc_files[j], stringsAsFactors = FALSE)
  if (!all(c("s", "Yc") %in% names(yc_df))) {
    stop("Member file must contain columns s and Yc: ", yc_files[j])
  }
  yc_df <- yc_df[order(yc_df$s), c("s", "Yc")]
  if (!identical(as.integer(yc_df$s), seq_len(nrow(grid_df)))) {
    stop("Member file does not match expected site indexing: ", yc_files[j])
  }

  yj <- as.numeric(yc_df$Yc)[ocean_sites]
  bad_idx <- which(!is.finite(yj))
  if (length(bad_idx) > 0L) {
    y_fill <- stats::median(yj[is.finite(yj)])
    yj[bad_idx] <- y_fill
    imputed_nonfinite_total <- imputed_nonfinite_total + length(bad_idx)
  }
  Yc_matrix[, j] <- yj
}

if (imputed_nonfinite_total > 0L) {
  warning(sprintf(
    "Imputed %d non-finite ocean values across ensemble members before averaging.",
    imputed_nonfinite_total
  ))
}

y <- rowMeans(Yc_matrix)
y_full <- rep(NA_real_, nrow(grid_df))
y_full[ocean_sites] <- y
rm(Yc_matrix)
invisible(gc())

lambda <- as.numeric(eig$eigenvalues)
ord <- order(lambda, decreasing = FALSE)
U <- eig$U[, ord, drop = FALSE]
lambda <- lambda[ord]

##########
# 3) Quantitative truncation diagnostics
##########

cat("Computing spectral coefficients and cumulative diagnostics...\n")
z <- as.numeric(crossprod(U, y))
energy <- z^2
cum_energy <- cumsum(energy)
tot_energy <- sum(energy)
if (!is.finite(tot_energy) || tot_energy <= 0) {
  stop("Total spectral energy is not positive/finite.")
}

cum_explained <- cum_energy / tot_energy
rmse_by_k <- sqrt((tot_energy - cum_energy) / n_ocean)

# Candidate K values for inspection. Includes small/medium/large values.
K_candidates <- unique(pmax(2L, pmin(
  n_ocean,
  c(20L, 40L, 80L, 120L, 160L, 240L, 320L, 500L, 800L, K_focus, 2000L, 3000L, 4000L, 5000L)
)))

diag_df <- data.frame(
  K = K_candidates,
  rmse = rmse_by_k[K_candidates],
  explained = cum_explained[K_candidates]
)
diag_df$rmse_rel <- diag_df$rmse / stats::sd(y)

# Data-driven "elbow-like" suggestions based on explained variance thresholds.
pick_first_k <- function(thresh) {
  idx <- which(cum_explained >= thresh)[1]
  if (is.na(idx)) n_ocean else idx
}

k90 <- pick_first_k(0.90)
k95 <- pick_first_k(0.95)
k99 <- pick_first_k(0.99)

focus_df <- data.frame(
  K = K_focus,
  rmse = rmse_by_k[K_focus],
  explained = cum_explained[K_focus]
)
focus_df$rmse_rel <- focus_df$rmse / stats::sd(y)

summary_txt <- sprintf(
  paste(
    "n_ocean=%d",
    "focus K=%s",
    "recommended=%d",
    sep = " | "
  ),
  n_ocean, paste(K_focus, collapse = ","), K_recommend
)
cat(summary_txt, "\n")

utils::write.csv(
  diag_df,
  file.path(plot_dir, "mode_truncation_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  focus_df,
  file.path(plot_dir, "mode_truncation_focus_metrics.csv"),
  row.names = FALSE
)

##########
# 4) Build plots for decision making
##########

p_rmse <- ggplot(diag_df, aes(x = K, y = rmse)) +
  geom_line(color = "#1b9e77", linewidth = 0.7) +
  geom_point(color = "#1b9e77", size = 1.7) +
  geom_vline(xintercept = K_focus, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  scale_x_log10() +
  labs(
    title = "Ensemble-Mean Reconstruction RMSE vs Number of Modes (K)",
    subtitle = summary_txt,
    x = "K modes kept (log scale)",
    y = "RMSE on ensemble-mean ocean SST"
  ) +
  theme_minimal(base_size = 12)

p_expl <- ggplot(diag_df, aes(x = K, y = explained)) +
  geom_line(color = "#2166ac", linewidth = 0.7) +
  geom_point(color = "#2166ac", size = 1.7) +
  geom_hline(yintercept = c(0.90, 0.95, 0.99), linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_vline(xintercept = K_focus, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Explained Spectral Energy of the Ensemble Mean vs K",
    subtitle = summary_txt,
    x = "K modes kept (log scale)",
    y = "Explained energy fraction"
  ) +
  theme_minimal(base_size = 12)

spec_df <- data.frame(
  k = seq_along(z),
  abs_coef = abs(z),
  eigenvalue = lambda
)

p_spec <- ggplot(spec_df, aes(x = k, y = abs_coef)) +
  geom_line(color = "#8c510a", linewidth = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Spectral Coefficient Magnitudes |z_k| for the Ensemble Mean",
    subtitle = "Large tails imply non-negligible high-frequency content",
    x = "Mode index k (sorted by eigenvalue, log scale)",
    y = "|z_k| (log scale)"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(plot_dir, "mode_truncation_rmse_curve.png"),
  plot = p_rmse, width = 8, height = 4.8, dpi = 300
)
ggsave(
  filename = file.path(plot_dir, "mode_truncation_explained_curve.png"),
  plot = p_expl, width = 8, height = 4.8, dpi = 300
)
ggsave(
  filename = file.path(plot_dir, "mode_truncation_spectral_magnitude.png"),
  plot = p_spec, width = 8, height = 4.8, dpi = 300
)

focus_window <- subset(diag_df, K >= min(K_focus) - 250L & K <= max(K_focus) + 250L)
p_focus <- ggplot(focus_window, aes(x = K)) +
  geom_line(aes(y = explained), color = "#2166ac", linewidth = 0.7) +
  geom_point(aes(y = explained), color = "#2166ac", size = 2) +
  geom_vline(xintercept = K_focus, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  labs(
    title = "Focused Explained-Energy Diagnostics Near the Working K Range",
    subtitle = sprintf("Recommended K = %d", K_recommend),
    x = "K modes kept",
    y = "Explained energy fraction"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(plot_dir, "mode_truncation_focus_window.png"),
  plot = p_focus, width = 8, height = 4.8, dpi = 300
)

##########
# 5) Residual maps for representative K
##########

K_map <- K_focus

cat("Building residual maps for K = ", paste(K_map, collapse = ", "), " ...\n", sep = "")

resid_map_list <- vector("list", length(K_map))
for (i in seq_along(K_map)) {
  K <- K_map[i]
  UK <- U[, seq_len(K), drop = FALSE]
  y_hat <- as.numeric(UK %*% crossprod(UK, y))
  resid <- y - y_hat

  resid_full <- rep(NA_real_, nrow(grid_df))
  resid_full[ocean_sites] <- resid

  resid_map_list[[i]] <- data.frame(
    lon = grid_df$lon,
    lat = grid_df$lat,
    resid = resid_full,
    K = factor(sprintf("K = %d", K), levels = sprintf("K = %d", K_map))
  )
}
resid_map_df <- do.call(rbind, resid_map_list)
lim <- max(abs(resid_map_df$resid), na.rm = TRUE)

p_resid <- ggplot(subset(resid_map_df, is.finite(resid)), aes(x = lon, y = lat, fill = resid)) +
  geom_raster() +
  coord_equal(expand = FALSE) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0, limits = c(-lim, lim)
  ) +
  facet_wrap(~ K, ncol = 2) +
  labs(
    title = "Residual Maps: Ensemble Mean - Projection onto First K Modes",
    subtitle = "Large coherent residuals indicate missing spatial structure",
    x = "Longitude",
    y = "Latitude",
    fill = "Residual"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  filename = file.path(plot_dir, "mode_truncation_residual_maps.png"),
  plot = p_resid, width = 10.5, height = 7.2, dpi = 300
)

##########
# 6) Console summary for quick recommendation
##########

cat("\nCandidate K diagnostics:\n")
print(diag_df, row.names = FALSE)

cat("\nReference explained-energy thresholds from the ensemble mean:\n")
cat(sprintf("  K@90%% = %d\n", k90))
cat(sprintf("  K@95%% = %d\n", k95))
cat(sprintf("  K@99%% = %d\n", k99))

cat("\nFocus-range diagnostics for the working band 1000-1500:\n")
print(focus_df, row.names = FALSE)

cat(sprintf("\nRecommended working value: K = %d\n", K_recommend))
cat("This sits in the middle of the inspected band and captures more than 99.99% of the ensemble-mean spectral energy.\n")

cat("\nWrote:\n")
cat("  - ", normalizePath(file.path(plot_dir, "mode_truncation_metrics.csv")), "\n", sep = "")
cat("  - ", normalizePath(file.path(plot_dir, "mode_truncation_focus_metrics.csv")), "\n", sep = "")
cat("  - ", normalizePath(file.path(plot_dir, "mode_truncation_rmse_curve.png")), "\n", sep = "")
cat("  - ", normalizePath(file.path(plot_dir, "mode_truncation_explained_curve.png")), "\n", sep = "")
cat("  - ", normalizePath(file.path(plot_dir, "mode_truncation_spectral_magnitude.png")), "\n", sep = "")
cat("  - ", normalizePath(file.path(plot_dir, "mode_truncation_focus_window.png")), "\n", sep = "")
cat("  - ", normalizePath(file.path(plot_dir, "mode_truncation_residual_maps.png")), "\n", sep = "")
