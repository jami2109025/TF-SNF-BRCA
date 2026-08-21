# ============================================================
# 07a_timing_probe.R
# Measure the actual per-operation cost on THIS machine, then
# project the bootstrap runtime for every candidate configuration.
# Runs in ~5 minutes.
# ============================================================
# WHY
# ------------------------------------------------------------
# Every runtime estimate so far has been extrapolated from the old
# pipeline's logs, and one of those extrapolations was already wrong:
# the distance cache was predicted at ~25 minutes and took 3.5. That
# error mattered in the dangerous direction -- if dist2 is cheap, then
# the 1150 core-seconds per replicate in the original 22.8-hour
# bootstrap were spent on SNF and spectralClustering instead, which
# would make an 8-view bootstrap far more expensive than projected.
#
# With a hard deadline there is no room to discover that at hour 40.
# This script measures the three operations that actually dominate --
# affinityMatrix, SNF, spectralClustering -- at full cohort size, and
# prints a runtime table for the real configuration choices.
#
# It also tests one legitimate lever: whether SNF's diffusion iteration
# count can be halved (t = 20 -> 10) without changing the clustering.
# SNF cost is linear in t, so if the partitions agree this halves the
# entire bootstrap for free. If they do not agree, the lever is not
# available and the script says so -- it is checked, not assumed.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SNFtool)
  library(cluster)
  library(mclust)
  library(survival)
  library(parallel)
})
source(here::here("R", "utils.R"))
source(here::here("R", "utils_benchmark.R"))

set_pipeline_seed()

SNF_K <- 20; SNF_SIGMA <- 0.5; SNF_T <- 20
n_cores <- max(1, parallel::detectCores() - 1)

cat("=========================================================\n")
cat("Timing probe -- measuring real per-operation cost\n")
cat("Cores available for the bootstrap:", n_cores, "\n")
cat("=========================================================\n\n")

# ------------------------------------------------------------------
# Load one expression view and one TF view at full size
# ------------------------------------------------------------------
vst_top2000 <- readRDS(here::here("data/processed/vst_top2000_genes.rds"))
tf_grid     <- readRDS(here::here("data/processed/tf_activity_grid.rds"))
tf_one      <- tf_grid[[if ("dorothea_AC__viper" %in% names(tf_grid))
                          "dorothea_AC__viper" else names(tf_grid)[1]]]

common <- intersect(colnames(vst_top2000), colnames(tf_one))
vst_top2000 <- vst_top2000[, common, drop = FALSE]
tf_one      <- tf_one[, common, drop = FALSE]
n <- length(common)
cat("Cohort size n =", n, "\n\n")

scale_view <- function(m) { s <- scale(t(m)); s[is.na(s)] <- 0; s }
expr_scaled <- scale_view(vst_top2000)
tf_scaled   <- scale_view(tf_one)

timeit <- function(expr, label) {
  t0 <- Sys.time()
  res <- force(expr)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  %-34s %8.1f s\n", label, el))
  list(result = res, secs = el)
}

cat("Measuring individual operations:\n")

t_dist_expr <- timeit(dist2_cached(expr_scaled), "dist2 (expression, 2000 feat)")
t_dist_tf   <- timeit(dist2_cached(tf_scaled),   "dist2 (TF view)")
D_expr <- t_dist_expr$result
D_tf   <- t_dist_tf$result

t_aff <- timeit(affinity_from_dist(D_expr, K = SNF_K, sigma = SNF_SIGMA), "affinityMatrix")
W_expr <- t_aff$result
W_tf   <- affinity_from_dist(D_tf, K = SNF_K, sigma = SNF_SIGMA)

t_snf <- timeit(SNFtool::SNF(list(W_expr, W_tf), K = SNF_K, t = SNF_T),
                paste0("SNF (t = ", SNF_T, ")"))
W_fused <- t_snf$result

t_spec <- timeit(SNFtool::spectralClustering(W_fused, K = 3), "spectralClustering")
cl3 <- t_spec$result

t_sil <- timeit(mean(cluster::silhouette(cl3, as.dist(1 - W_fused))[, "sil_width"]),
                "silhouette (native graph)")

# ------------------------------------------------------------------
# Is the t = 20 -> 10 lever available?
# ------------------------------------------------------------------
cat("\nTesting whether SNF diffusion iterations can be halved:\n")
t_snf10 <- timeit(SNFtool::SNF(list(W_expr, W_tf), K = SNF_K, t = 10), "SNF (t = 10)")
cl3_t10 <- SNFtool::spectralClustering(t_snf10$result, K = 3)
cl2_t10 <- SNFtool::spectralClustering(t_snf10$result, K = 2)
cl2     <- SNFtool::spectralClustering(W_fused, K = 2)

ari_t10_k2 <- mclust::adjustedRandIndex(cl2, cl2_t10)
ari_t10_k3 <- mclust::adjustedRandIndex(cl3, cl3_t10)
cat(sprintf("  ARI(t=20 vs t=10):  k=2: %.4f   k=3: %.4f\n", ari_t10_k2, ari_t10_k3))

t_lever <- (ari_t10_k2 > 0.99 && ari_t10_k3 > 0.99)
if (t_lever) {
  cat("  -> LEVER AVAILABLE. t=10 reproduces the t=20 partition (ARI > 0.99).\n")
  cat("     Halving SNF_T halves the bootstrap. Report the change and the ARI\n")
  cat("     check in the methods section; do not make it silently.\n")
} else {
  cat("  -> LEVER NOT AVAILABLE. t=10 changes the partition. Keep SNF_T = 20.\n")
}

# ------------------------------------------------------------------
# Cost model
# ------------------------------------------------------------------
snf_s  <- t_snf$secs
aff_s  <- t_aff$secs
spec_s <- t_spec$secs
sil_s  <- t_sil$secs

# Per view per replicate inside the bootstrap, with the corrected loop
# order (network built ONCE, then clustered at each k):
#   1 affinity (TF) + 1 SNF + length(K_VALUES) spectralClustering
# Endpoint evaluation (ARI/NMI/Cox) is cheap relative to these and is
# folded in as a 10% allowance.
per_view_rep_s <- (aff_s + snf_s + 2 * spec_s) * 1.10

# Per view in the point-estimate stage:
#   2 affinity + 1 SNF + 6 spectralClustering + 6 silhouette
per_view_point_s <- (2 * aff_s + snf_s + 6 * spec_s + 6 * sil_s) * 1.10

cat("\n=========================================================\n")
cat("COST MODEL\n")
cat("=========================================================\n")
cat(sprintf("  per view, per bootstrap replicate : %6.1f core-seconds\n", per_view_rep_s))
cat(sprintf("  per view, point-estimate stage    : %6.1f core-seconds\n", per_view_point_s))

boot_hours <- function(n_views, n_boot, t_halved = FALSE) {
  cost <- if (t_halved) (aff_s + snf_s / 2 + 2 * spec_s) * 1.10 else per_view_rep_s
  n_views * n_boot * cost / n_cores / 3600
}

cat("\nPROJECTED BOOTSTRAP RUNTIME (hours)\n\n")
grid <- expand.grid(n_views = c(2, 4, 6, 8), n_boot = c(200, 300, 400, 500))
proj <- grid %>%
  mutate(hours_t20 = mapply(function(v, b) boot_hours(v, b, FALSE), n_views, n_boot),
         hours_t10 = mapply(function(v, b) boot_hours(v, b, TRUE),  n_views, n_boot)) %>%
  arrange(n_views, n_boot)

print(proj %>% mutate(across(starts_with("hours"), ~ round(.x, 1))), n = 100)

BUDGET_HOURS <- 40
cat("\nConfigurations fitting a", BUDGET_HOURS, "hour budget",
    if (t_lever) "(using t = 10, the verified lever):\n" else "(t = 20):\n")
feasible <- proj %>%
  mutate(h = if (t_lever) hours_t10 else hours_t20) %>%
  filter(h <= BUDGET_HOURS) %>%
  arrange(desc(n_views * n_boot))
if (nrow(feasible) == 0) {
  cat("  NONE. Reduce the grid further, or accept point estimates without CIs\n")
  cat("  for some views and say so explicitly in the manuscript.\n")
} else {
  print(feasible %>% mutate(h = round(h, 1)) %>% dplyr::select(n_views, n_boot, h), n = 20)
  best <- feasible[1, ]
  cat(sprintf("\n  RECOMMENDED: VIEWS_FOR_BOOTSTRAP = %d views, N_BOOT = %d  (~%.1f h)\n",
              best$n_views, best$n_boot, best$h))
  if (t_lever) cat("               and set SNF_T <- 10 in 07_benchmark_models.R\n")
}

cat("\nAlso note the point-estimate stage is SERIAL in the current script:\n")
cat(sprintf("  %d views would take ~%.1f hours; %d views ~%.1f hours.\n",
            6, 6 * per_view_point_s / 3600, 30, 30 * per_view_point_s / 3600))
cat("  Parallelise it (see the patched 07) if the full grid is large.\n")
cat("=========================================================\n")
