# ============================================================
# 06_snf_clustering.R  (publication-ready revision, FINAL)
# Similarity Network Fusion: expression view + TF-activity view
# ============================================================
# CHANGES in this revision (post field-testing on the real cohort):
#   - Consensus PAC and the hyperparameter sweep are now PARALLELIZED
#     across cores (parallel::makeCluster/parLapply), since the
#     sequential version took >25 hours on real data. Seeding is
#     preserved per-repetition (PIPELINE_SEED + offset + b), so results
#     are bit-identical to a sequential run -- only wall-clock time
#     changes. Set N_CORES_OVERRIDE below if you want to control this
#     explicitly (e.g. on a shared machine).
#   - Fixed a bug where `sweep_grid` was not included in the sweep's
#     clusterExport() varlist, which would have crashed the
#     hyperparameter sweep AFTER the multi-hour consensus PAC step
#     completed.
#   - Fused-network heatmaps are now colored using data-driven
#     percentile breaks (not a fixed 0-0.5 scale), since the real
#     affinity values in this cohort top out far below 0.5 (~0.06),
#     which made the original fixed-scale heatmap appear uniformly
#     flat and hid genuine block-diagonal cluster structure.
#   - NEW: automated batch-driver sensitivity check. If the formal
#     cluster-vs-batch association test (below) is significant at
#     FDR<0.05 for possible_plate, the script now automatically
#     identifies the single plate contributing the largest chi-square
#     standardized residual, re-clusters the fused network EXCLUDING
#     that plate's samples, and reports the ARI between the original
#     and plate-excluded cluster assignments. This distinguishes
#     "associated with a batch proxy" from "driven by a batch proxy"
#     and is saved as a tracked output rather than only having existed
#     as ad hoc console commands.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SNFtool)
  library(cluster)
  library(pheatmap)
  library(parallel)
})
source(here::here("R", "utils.R"))

script_name <- "06_snf_clustering"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

# Set this to a fixed number if you want to control core usage
# explicitly (e.g. on a shared machine); NULL = auto-detect - 1.
N_CORES_OVERRIDE <- NULL
n_cores <- if (is.null(N_CORES_OVERRIDE)) max(1, parallel::detectCores() - 1) else N_CORES_OVERRIDE
cat("Using", n_cores, "cores for parallelized steps in this script.\n")

# -----------------------------
# Load and align inputs (expression view + TF-activity view)
# -----------------------------
vst_top2000     <- readRDS(here::here("data/processed/vst_top2000_genes.rds"))
tf_activity     <- readRDS(here::here("data/processed/tf_activity_viper_AC_primary.rds"))
sample_metadata <- readRDS(here::here("data/processed/sample_metadata_matched.rds"))

aligned <- align_samples(expr = vst_top2000, tf = tf_activity, meta = sample_metadata)
vst_top2000 <- aligned$expr
tf_activity <- aligned$tf
sample_metadata <- aligned$meta
stopifnot(identical(colnames(vst_top2000), colnames(tf_activity)))
stopifnot(identical(colnames(vst_top2000), rownames(sample_metadata)))

cat("Feature-space dimensionality by view -- expression:", nrow(vst_top2000),
    "genes; TF activity:", nrow(tf_activity), "TFs. (Reported explicitly since",
    "the two views are asymmetric in dimensionality; SNF fuses PER-VIEW\n",
    "sample-similarity graphs, so this does not directly bias fusion the way\n",
    "naive feature concatenation would, but is disclosed for transparency.)\n")

# -----------------------------
# Per-feature z-scoring of each view, then NA -> 0 imputation post-scaling
# -----------------------------
expr_scaled <- scale(t(vst_top2000))
tf_scaled   <- scale(t(tf_activity))
expr_scaled[is.na(expr_scaled)] <- 0
tf_scaled[is.na(tf_scaled)]     <- 0

cat("Expression view:", dim(expr_scaled), "  TF-activity view:", dim(tf_scaled), "\n")

# -----------------------------
# Core SNF construction, reused for the base run, the subsampling-based
# consensus check, AND the hyperparameter sensitivity sweep.
# -----------------------------
build_fused_network <- function(expr_mat, tf_mat, K = 20, sigma = 0.5, t_iter = 20) {
  D_expr <- SNFtool::dist2(as.matrix(expr_mat), as.matrix(expr_mat))
  D_tf   <- SNFtool::dist2(as.matrix(tf_mat),   as.matrix(tf_mat))
  W_expr <- SNFtool::affinityMatrix(D_expr, K = K, sigma = sigma)
  W_tf   <- SNFtool::affinityMatrix(D_tf,   K = K, sigma = sigma)
  W_fused <- SNFtool::SNF(list(W_expr, W_tf), K = K, t = t_iter)
  list(W_expr = W_expr, W_tf = W_tf, W_fused = W_fused)
}

BASE_K <- 20; BASE_SIGMA <- 0.5; BASE_T <- 20
base_net <- build_fused_network(expr_scaled, tf_scaled, K = BASE_K, sigma = BASE_SIGMA, t_iter = BASE_T)
W_fused  <- base_net$W_fused
rownames(W_fused) <- colnames(W_fused) <- colnames(vst_top2000)

# -----------------------------
# Eigengap (k-selection diagnostic)
# -----------------------------
calculate_eigengap <- function(W, max_k = 6) {
  n <- nrow(W); d <- rowSums(W)
  D_inv_sqrt <- diag(1 / sqrt(d + 1e-10))
  L_norm <- diag(n) - D_inv_sqrt %*% W %*% D_inv_sqrt
  eigen_values <- sort(eigen(L_norm, symmetric = TRUE, only.values = TRUE)$values)
  gaps <- eigen_values[2:(max_k + 1)] - eigen_values[1:max_k]
  tibble(
    k = 1:max_k,
    eigenvalue_k = eigen_values[1:max_k],
    eigengap_after_k = gaps,
    interpretation = paste0("Gap between eigenvalue ", 1:max_k, " and ", 2:(max_k + 1),
                            "; larger values support k=", 1:max_k, " clusters")
  )
}
eigengap_table <- calculate_eigengap(W_fused, max_k = 6)
cat("\nEigengap table:\n"); print(eigengap_table %>% dplyr::select(-interpretation))

# -----------------------------
# Silhouette on the fused graph (graph-distance proxy)
# -----------------------------
calculate_network_silhouette <- function(W, k) {
  clusters <- SNFtool::spectralClustering(W, K = k)
  distance_matrix <- as.dist(1 - W)
  sil <- cluster::silhouette(clusters, distance_matrix)
  mean(sil[, "sil_width"])
}
k_values <- 2:6
silhouette_table <- tibble(
  k = k_values,
  network_silhouette = sapply(k_values, function(k) calculate_network_silhouette(W_fused, k))
)
cat("\nNetwork silhouette by k (graph-distance proxy, NOT feature-space silhouette):\n")
print(silhouette_table)

# -----------------------------
# Consensus PAC via SUBSAMPLING (Monti et al. 2003 convention), PARALLELIZED
# -----------------------------
consensus_pac_subsampling <- function(expr_mat, tf_mat, k, n_reps = 50, frac = 0.8,
                                      K = BASE_K, sigma = BASE_SIGMA, t_iter = BASE_T) {
  n <- nrow(expr_mat); n_sub <- round(frac * n)
  
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterEvalQ(cl, library(SNFtool))
  parallel::clusterExport(
    cl,
    varlist = c("expr_mat", "tf_mat", "n", "n_sub", "K", "sigma", "t_iter",
                "k", "build_fused_network", "PIPELINE_SEED"),
    envir = environment()
  )
  
  rep_results <- parallel::parLapply(cl, seq_len(n_reps), function(b) {
    set.seed(PIPELINE_SEED + 1000 + b)
    idx <- sample(seq_len(n), n_sub)
    sub_ids <- rownames(expr_mat)[idx]
    net_b <- build_fused_network(expr_mat[idx, ], tf_mat[idx, ], K = K, sigma = sigma, t_iter = t_iter)
    cl_b  <- SNFtool::spectralClustering(net_b$W_fused, K = k)
    names(cl_b) <- sub_ids
    list(sub_ids = sub_ids, cl_b = cl_b)
  })
  
  consensus_mat <- matrix(0, n, n, dimnames = list(rownames(expr_mat), rownames(expr_mat)))
  count_mat <- matrix(0, n, n, dimnames = dimnames(consensus_mat))
  for (r in rep_results) {
    same_cluster <- outer(r$cl_b, r$cl_b, FUN = "==")
    consensus_mat[r$sub_ids, r$sub_ids] <- consensus_mat[r$sub_ids, r$sub_ids] + same_cluster
    count_mat[r$sub_ids, r$sub_ids]     <- count_mat[r$sub_ids, r$sub_ids] + 1
  }
  consensus_frac <- consensus_mat / pmax(count_mat, 1)
  diag(consensus_frac) <- NA
  vals <- consensus_frac[upper.tri(consensus_frac)]
  vals <- vals[!is.na(vals)]
  mean(vals > 0.1 & vals < 0.9)
}

cat("\nRunning FULL-PIPELINE, subsampling-based consensus PAC (parallelized,",
    n_cores, "cores) -- may take a while.\n")

pac_table <- tibble(
  k = k_values,
  consensus_PAC = sapply(k_values, function(k) {
    consensus_pac_subsampling(expr_scaled, tf_scaled, k = k, n_reps = 50)
  })
)
cat("\nConsensus PAC by k (subsampling-based; LOWER = more stable):\n")
print(pac_table)

k_selection_metrics <- eigengap_table %>%
  dplyr::select(k, eigenvalue_k, eigengap_after_k) %>%
  inner_join(silhouette_table, by = "k") %>%
  inner_join(pac_table, by = "k")
write_csv(k_selection_metrics, here::here("results/tables/table_k_selection_metrics.csv"))
cat("\nFull k-selection metrics table:\n"); print(k_selection_metrics)

# -----------------------------
# Hyperparameter sensitivity sweep (K, sigma), PARALLELIZED
# (FIX applied: sweep_grid is now correctly included in clusterExport)
# -----------------------------
cat("\nRunning SNF hyperparameter sensitivity sweep (K x sigma), parallelized...\n")

sweep_grid <- expand.grid(K_param = c(10, 20, 30), sigma_param = c(0.3, 0.5, 0.8))

base_clusters_by_k <- lapply(c(2, 3), function(k) SNFtool::spectralClustering(W_fused, K = k))
names(base_clusters_by_k) <- c("k2", "k3")

cl_sweep <- parallel::makeCluster(n_cores)
parallel::clusterEvalQ(cl_sweep, library(SNFtool))
parallel::clusterExport(cl_sweep, varlist = c("expr_scaled", "tf_scaled", "BASE_T", "sweep_grid",
                                              "build_fused_network", "base_clusters_by_k"))
sweep_results_list <- parallel::parLapply(cl_sweep, seq_len(nrow(sweep_grid)), function(i) {
  K_param <- sweep_grid$K_param[i]; sigma_param <- sweep_grid$sigma_param[i]
  net_sw <- build_fused_network(expr_scaled, tf_scaled, K = K_param, sigma = sigma_param, t_iter = BASE_T)
  tibble::tibble(
    K_param = K_param, sigma_param = sigma_param,
    ari_k2 = mclust::adjustedRandIndex(SNFtool::spectralClustering(net_sw$W_fused, K = 2), base_clusters_by_k$k2),
    ari_k3 = mclust::adjustedRandIndex(SNFtool::spectralClustering(net_sw$W_fused, K = 3), base_clusters_by_k$k3)
  )
})
parallel::stopCluster(cl_sweep)
sensitivity_results <- dplyr::bind_rows(sweep_results_list)

cat("\nHyperparameter sensitivity (ARI vs. base K=20, sigma=0.5 clustering):\n")
print(sensitivity_results)
write_csv(sensitivity_results, here::here("results/tables/table_snf_hyperparameter_sensitivity.csv"))

# -----------------------------
# Final cluster assignments -- BOTH k=2 and k=3, co-equal
# -----------------------------
final_clusters_k2 <- SNFtool::spectralClustering(W_fused, K = 2)
final_clusters_k3 <- SNFtool::spectralClustering(W_fused, K = 3)

cluster_labels_k2 <- tibble(
  full_barcode = colnames(vst_top2000), k = 2L,
  cluster_int = final_clusters_k2, SNF_cluster = paste0("SNF_C", final_clusters_k2)
)
cluster_labels_k3 <- tibble(
  full_barcode = colnames(vst_top2000), k = 3L,
  cluster_int = final_clusters_k3, SNF_cluster = paste0("SNF_C", final_clusters_k3)
)

cat("\nk=2 cluster sizes:\n"); print(table(cluster_labels_k2$SNF_cluster))
cat("\nk=3 cluster sizes:\n"); print(table(cluster_labels_k3$SNF_cluster))

ari_k2_vs_k3 <- mclust::adjustedRandIndex(final_clusters_k2, final_clusters_k3)
cat("\nARI between k=2 and k=3 solutions:", round(ari_k2_vs_k3, 3), "\n")

# -----------------------------
# Formal cluster-vs-batch association test
# -----------------------------
batch_available <- all(c("possible_plate", "possible_center") %in% colnames(sample_metadata))
batch_association_table <- tibble()

if (batch_available) {
  batch_test_k2_center <- test_cluster_batch_association(
    cluster_labels_k2$SNF_cluster, sample_metadata$possible_center, seed_offset = 1400)
  batch_test_k3_center <- test_cluster_batch_association(
    cluster_labels_k3$SNF_cluster, sample_metadata$possible_center, seed_offset = 1401)
  batch_test_k2_plate <- test_cluster_batch_association(
    cluster_labels_k2$SNF_cluster, sample_metadata$possible_plate, seed_offset = 1402)
  batch_test_k3_plate <- test_cluster_batch_association(
    cluster_labels_k3$SNF_cluster, sample_metadata$possible_plate, seed_offset = 1403)
  
  batch_association_table <- tibble(
    k = c(2, 3, 2, 3),
    batch_variable = c("possible_center", "possible_center", "possible_plate", "possible_plate"),
    test = c(batch_test_k2_center$test, batch_test_k3_center$test,
             batch_test_k2_plate$test, batch_test_k3_plate$test),
    p_value = c(batch_test_k2_center$p_value, batch_test_k3_center$p_value,
                batch_test_k2_plate$p_value, batch_test_k3_plate$p_value)
  ) %>%
    mutate(FDR = p.adjust(p_value, method = "BH"))
  
  cat("\n\nCluster-vs-technical-batch association:\n")
  print(batch_association_table)
  write_csv(batch_association_table, here::here("results/tables/table_cluster_batch_association.csv"))
  
  if (any(batch_association_table$FDR < 0.05, na.rm = TRUE)) {
    warning(
      "SNF cluster assignment shows a statistically significant association ",
      "with a technical batch proxy at FDR<0.05. See ",
      "table_batch_driver_sensitivity_check.csv for an automated exclusion ",
      "sensitivity test identifying whether this is DRIVEN by, or merely ",
      "ASSOCIATED with, the top-residual batch category."
    )
  }
} else {
  warning("possible_plate/possible_center not found in sample_metadata; ",
          "run the updated 04_quality_control.R first. Batch association NOT tested.")
}

# -----------------------------
# NEW: automated batch-driver sensitivity check.
# If any k x batch_variable combination is significant at FDR<0.05,
# identify the single category (e.g. one plate) contributing the
# largest chi-square standardized residual, exclude its samples, and
# re-cluster to test whether the ORIGINAL partition is reproduced.
# High ARI => associated but not driven by that category (as found for
# A00Z during development on the TCGA-BRCA cohort: ARI=0.979 for both
# k=2 and k=3). Low ARI => the category is materially shaping the
# partition and warrants deeper investigation / discussion as a
# limitation.
# -----------------------------
find_top_residual_category <- function(cluster_vector, batch_vector) {
  tab <- table(cluster_vector, batch_vector)
  tab <- tab[, colSums(tab) > 0, drop = FALSE]
  chi <- suppressWarnings(chisq.test(tab))
  resid <- chi$stdres
  idx <- which(abs(resid) == max(abs(resid)), arr.ind = TRUE)[1, ]
  list(category = colnames(tab)[idx["batch_vector"]], std_resid = resid[idx["cluster_vector"], idx["batch_vector"]])
}

batch_sensitivity_results <- tibble()
if (batch_available && nrow(batch_association_table) > 0 &&
    any(batch_association_table$FDR < 0.05, na.rm = TRUE)) {
  
  sig_rows <- batch_association_table %>% filter(FDR < 0.05)
  
  batch_sensitivity_results <- purrr::pmap_dfr(sig_rows, function(k, batch_variable, test, p_value, FDR) {
    cluster_labels_this_k <- if (k == 2) cluster_labels_k2 else cluster_labels_k3
    batch_vec <- sample_metadata[[batch_variable]]
    
    top_cat <- find_top_residual_category(cluster_labels_this_k$SNF_cluster, batch_vec)
    cat("\nTop residual category for k=", k, ", ", batch_variable, ": '", top_cat$category,
        "' (std. residual = ", round(top_cat$std_resid, 2), ")\n", sep = "")
    
    excl_samples <- rownames(sample_metadata)[batch_vec == top_cat$category & !is.na(batch_vec)]
    keep_samples <- setdiff(colnames(W_fused), excl_samples)
    
    W_sub <- W_fused[keep_samples, keep_samples]
    cl_sub <- SNFtool::spectralClustering(W_sub, K = k)
    names(cl_sub) <- keep_samples
    
    orig_labels <- cluster_labels_this_k %>%
      filter(full_barcode %in% keep_samples) %>%
      arrange(match(full_barcode, keep_samples)) %>%
      pull(cluster_int)
    
    ari_sensitivity <- mclust::adjustedRandIndex(cl_sub, orig_labels)
    
    tibble(
      k = k, batch_variable = batch_variable, top_category = top_cat$category,
      std_residual = top_cat$std_resid, n_excluded = length(excl_samples),
      n_remaining = length(keep_samples), ari_original_vs_excluded = ari_sensitivity,
      interpretation = ifelse(
        ari_sensitivity >= 0.8,
        "High ARI: cluster structure reproduced without this category -- ASSOCIATED but not DRIVEN by it.",
        "Lower ARI: cluster structure changes materially without this category -- requires further investigation."
      )
    )
  })
  
  cat("\n\n=== Automated batch-driver sensitivity check ===\n")
  print(batch_sensitivity_results)
  write_csv(batch_sensitivity_results, here::here("results/tables/table_batch_driver_sensitivity_check.csv"))
} else {
  cat("\nNo FDR-significant batch association found; batch-driver sensitivity check skipped.\n")
}

# Combined object consumed by every downstream script
final_snf_cluster_labels <- list(
  k2 = cluster_labels_k2,
  k3 = cluster_labels_k3,
  k_selection_metrics = k_selection_metrics,
  hyperparameter_sensitivity = sensitivity_results,
  batch_association = batch_association_table,
  batch_driver_sensitivity = batch_sensitivity_results,
  ari_k2_vs_k3 = ari_k2_vs_k3,
  base_params = list(K = BASE_K, sigma = BASE_SIGMA, t_iter = BASE_T),
  note = paste(
    "Both k=2 and k=3 SNF clusterings are reported as CO-EQUAL, parallel",
    "solutions throughout this pipeline. See table_k_selection_metrics.csv,",
    "table_cluster_batch_association.csv, and",
    "table_batch_driver_sensitivity_check.csv for full diagnostics."
  )
)
saveRDS(final_snf_cluster_labels, here::here("results/objects/final_snf_cluster_labels.rds"))
saveRDS(W_fused, here::here("results/objects/W_fused.rds"))

write_csv(bind_rows(cluster_labels_k2, cluster_labels_k3),
          here::here("results/tables/table_snf_cluster_labels_k2_k3.csv"))

# -----------------------------
# Heatmap of the fused network, ordered by cluster, with DATA-DRIVEN
# color breaks (99th percentile of off-diagonal values), since a fixed
# 0-0.5 scale hides real structure when true values top out near 0.06.
# -----------------------------
val_99 <- quantile(W_fused[upper.tri(W_fused)], 0.99)
color_breaks <- seq(0, val_99, length.out = 101)
cat("\nHeatmap color scale capped at 99th percentile of off-diagonal values:",
    round(val_99, 5), "\n")

for (kk in list(list(lab = cluster_labels_k2, name = "k2"), list(lab = cluster_labels_k3, name = "k3"))) {
  ord <- order(kk$lab$SNF_cluster)
  ann <- data.frame(SNF_cluster = kk$lab$SNF_cluster[ord], row.names = kk$lab$full_barcode[ord])
  png(here::here("results/figures", paste0("figure_snf_fused_heatmap_", kk$name, ".png")),
      width = 1400, height = 1200, res = 150)
  pheatmap(
    W_fused[kk$lab$full_barcode[ord], kk$lab$full_barcode[ord]],
    cluster_rows = FALSE, cluster_cols = FALSE,
    annotation_row = ann, annotation_col = ann,
    show_rownames = FALSE, show_colnames = FALSE,
    breaks = color_breaks,
    main = paste0("Fused SNF similarity network, ordered by ", kk$name, " clusters")
  )
  dev.off()
}

log_session_info(script_name, key_packages = c("SNFtool", "mclust", "cluster", "parallel"))
cat("\n\u2713 06_snf_clustering.R complete. Both k=2 and k=3 solutions saved;",
    "batch association and driver-sensitivity checks complete.\n")

close_logger(log_con, script_name)
