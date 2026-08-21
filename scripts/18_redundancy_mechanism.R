# ============================================================
# 18_redundancy_mechanism.R   (NEW -- explanatory, not confirmatory)
#
# WHY THIS SCRIPT EXISTS
# ------------------------------------------------------------
# Scripts 07 and 14 established WHAT happens: adding a TF-activity
# view does not improve agreement with PAM50 (inconclusive, the
# discovery cohort is underpowered) and does not improve adjusted
# prognostic discrimination in METABRIC (equivalent to zero, 29/29
# views, both k, with MDD four- to seventeenfold below the SESOI).
#
# Neither script explains WHY. A referee -- and a supervisor asking
# "what is your contribution?" -- will reasonably ask whether the null
# is a property of TF biology, of the SNF machinery, or of the
# experimental design. A negative result that offers a MECHANISM is
# citable; one that only reports a measurement is not.
#
# The hypothesis tested here is redundancy:
#
#   TF activity is a weighted average of the expression of the genes a
#   TF regulates. Those target genes are drawn from the same
#   transcriptome that supplies the expression view. If the top-2000
#   variance-selected genes already span the space in which TF activity
#   lives, then the TF view is a near-linear re-encoding of information
#   the model already has -- and the observed null is not a surprise
#   about regulatory biology, it is a predictable consequence of the
#   feature space.
#
# Three levels of evidence, each answering a different objection:
#
#   PART A -- FEATURE LEVEL. Cross-validated R^2 for predicting each
#     TF's activity from principal components of the expression view.
#     Out-of-sample, so it cannot be inflated by overfitting; a
#     sample-shuffled negative control confirms the estimator returns
#     ~0 when there is nothing to predict.
#     Objection answered: "maybe the views are unrelated and SNF just
#     failed to use the TF one."
#
#   PART B -- GEOMETRY LEVEL. SNF does not see features; it sees a
#     similarity graph. So the feature-level result is not sufficient.
#     Here we measure how far fusion actually MOVES the geometry:
#     correlation between the expression and TF distance matrices,
#     between the fused and expression-only affinity matrices, and --
#     most interpretably -- the ARI between the fused partition and the
#     expression-only partition.
#     Objection answered: "R^2 on features says nothing about what a
#     graph-fusion method does."
#
#   PART C -- MACHINERY LEVEL. Does the FUSION algorithm add anything
#     over simply gluing the feature matrices together? This is a
#     separate claim from the added-value question -- it is about SNF,
#     not about TF activity -- and the data to answer it are already in
#     table_benchmark_grid_k2_k3.csv from script 07. No recomputation.
#     Objection answered: "your null is about SNF, not about the view."
#
# WHAT THIS SCRIPT IS NOT
# ------------------------------------------------------------
# Every quantity here is DESCRIPTIVE and POST HOC. Nothing is tested
# against a pre-specified SESOI, because none was registered for these
# comparisons, and inventing one now -- after seeing scripts 07 and 14
# -- would be exactly the forking-path error the whole framework exists
# to prevent. These numbers explain the pre-specified results; they do
# not add new confirmatory claims, and the output labels them as such.
#
# COST -- READ THIS BEFORE LAUNCHING
# Single-threaded by design. No PSOCK cluster, so none of the failure
# modes that cost script 07 seven hours can occur. Peak memory is one
# 1095 x 1095 double (~9.6 MB) plus the expression view; views are
# processed one at a time and freed.
#
# Part A is fast (~2-5 min for all 29 views: one QR per fold serves
# every TF at once).
#
# Part B is split. B1 (all 29 views, no SNF) is ~3 minutes total.
# B2 needs a real SNF run per view -- 210 GFLOP at n=1095, t=20 --
# which script 07's own log prices at 15-20 minutes per view on this
# machine. B2 therefore runs on 5 pre-specified views only: ~1.5 hours.
#
# Running B2 on all 29 would cost ~8 hours to restate what the first
# five already show, and nothing in the thesis depends on it.
# RUN_ALL_VIEWS_IN_B2 exists for completeness and should stay FALSE.
#
# Total expected: roughly 1.5-2 hours, checkpointed throughout.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SNFtool)
})
source(here::here("R", "utils.R"))
source(here::here("R", "utils_benchmark.R"))

script_name <- "18_redundancy_mechanism"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

# ------------------------------------------------------------------
# Configuration. Fixed by rule, not tuned -- see note on PC selection.
# ------------------------------------------------------------------
SNF_K        <- 20
SNF_SIGMA    <- 0.5
SNF_T        <- 20
K_VALUES     <- c(2, 3)
N_FOLDS      <- 5
PC_VAR_TARGET <- 0.80   # retain PCs explaining this share of expression variance
MAX_PCS       <- 200    # hard cap, keeps n/p comfortable for CV regression
PC_SENSITIVITY <- 50    # second, much smaller PC count -- see note at Part A

# PART B PRIORITY SET. These are the views the thesis needs: the two
# carried through script 07's bootstrap, their degree-matched permuted
# counterparts (so the geometry comparison has its own control), and one
# extra regulon source for breadth. They are processed first so that an
# interrupted run still yields a complete, reportable Part B.
GEOMETRY_PRIORITY <- c(
  "dorothea_AC__viper", "dorothea_AC_edgeperm__viper",
  "collectri__viper",   "collectri_edgeperm__viper",
  "dorothea_AB__viper"
)
# Leave FALSE. TRUE adds ~7 hours to restate a result the priority set
# already shows, and nothing in the thesis depends on it.
RUN_ALL_VIEWS_IN_B2 <- FALSE

cat("\n============================================================\n")
cat("18_redundancy_mechanism -- MECHANISTIC / DESCRIPTIVE ANALYSIS\n")
cat("All quantities below are post hoc and explanatory. No\n")
cat("confirmatory claim is made from this script.\n")
cat("============================================================\n\n")

set_blas_threads_single()

# ------------------------------------------------------------------
# Load and align, exactly as script 07 does, so the geometry examined
# here is the SAME geometry the benchmark used. Any divergence in
# sample set or scaling would make Part B measure a different object
# from the one the added-value estimates came from.
# ------------------------------------------------------------------
vst_top2000     <- readRDS(here::here("data/processed/vst_top2000_genes.rds"))
sample_metadata <- readRDS(here::here("data/processed/sample_metadata_matched.rds"))
tf_grid         <- readRDS(here::here("data/processed/tf_activity_grid.rds"))

common_samples <- Reduce(intersect, c(
  list(colnames(vst_top2000), rownames(sample_metadata)),
  lapply(tf_grid, colnames)
))
stopifnot(length(common_samples) > 100)
cat("Canonical sample set:", length(common_samples), "patients\n")
cat("TF-activity grid:", length(tf_grid), "views\n\n")

vst_top2000 <- vst_top2000[, common_samples, drop = FALSE]
tf_grid     <- lapply(tf_grid, function(m) m[, common_samples, drop = FALSE])

scale_view <- function(mat) {
  s <- scale(t(mat))          # samples in rows
  s[is.na(s)] <- 0
  s
}
expr_scaled <- scale_view(vst_top2000)
n_samples <- nrow(expr_scaled)

# ==================================================================
# PART A -- FEATURE-LEVEL REDUNDANCY
# ==================================================================
cat("=== PART A: how much of TF activity is recoverable from the\n")
cat("    expression features already in the model? ===\n\n")

# PC count is chosen by a fixed variance rule rather than by tuning.
# This matters: if the number of PCs were selected to maximise R^2, the
# headline "TF activity is X% predictable" would be a fitted quantity
# rather than a measured one. The rule is stated, applied once, and the
# resulting count is reported so a reader can judge it.
pca <- prcomp(expr_scaled, center = FALSE, scale. = FALSE)
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
n_pc <- min(which(cumsum(var_explained) >= PC_VAR_TARGET)[1], MAX_PCS)
if (is.na(n_pc)) n_pc <- MAX_PCS
X_pc <- pca$x[, seq_len(n_pc), drop = FALSE]

cat("Expression PCs retained:", n_pc,
    sprintf("(%.1f%% of variance; rule = first PCs reaching %.0f%%, cap %d)\n",
            100 * cumsum(var_explained)[n_pc], 100 * PC_VAR_TARGET, MAX_PCS))
cat("Samples:", n_samples, " -> samples per predictor:",
    round(n_samples / n_pc, 1), "\n\n")

#' Out-of-sample R^2 by K-fold cross-validation.
#'
#' In-sample R^2 with n_pc predictors would be inflated by roughly
#' n_pc/n even for pure noise (here ~ 0.05-0.2), which is the same
#' order as an interesting effect. Cross-validation removes that
#' inflation entirely: predictions for each fold come from a model that
#' never saw it. A value near 0 therefore means "not predictable", and
#' negative values are possible and meaningful (worse than the mean).
#' Y is a MATRIX (samples x TFs). All TFs share one design matrix, so a
#' single QR decomposition per fold serves every TF at once. Fitting them
#' one at a time would repeat the same decomposition ~250 times per view
#' and turn a two-minute job into a half-hour one for no gain.
cv_r2_matrix <- function(Y, X, folds) {
  Y <- as.matrix(Y)
  pred <- matrix(NA_real_, nrow(Y), ncol(Y))
  for (f in unique(folds)) {
    tr <- folds != f; te <- !tr
    fit <- tryCatch(
      stats::lm.fit(x = cbind(1, X[tr, , drop = FALSE]), y = Y[tr, , drop = FALSE]),
      error = function(e) NULL)
    if (is.null(fit)) return(rep(NA_real_, ncol(Y)))
    cf <- fit$coefficients
    cf[is.na(cf)] <- 0
    pred[te, ] <- cbind(1, X[te, , drop = FALSE]) %*% cf
  }
  vapply(seq_len(ncol(Y)), function(j) {
    y <- Y[, j]; p <- pred[, j]
    ok <- is.finite(p) & is.finite(y)
    if (sum(ok) < 10 || sd(y[ok]) < 1e-10) return(NA_real_)
    1 - sum((y[ok] - p[ok])^2) / sum((y[ok] - mean(y[ok]))^2)
  }, numeric(1))
}

set.seed(20260815)
folds <- sample(rep_len(seq_len(N_FOLDS), n_samples))

# Only the REAL regulon views are informative for the redundancy claim;
# edge-permuted views are carried through as a comparison because a
# randomised regulon is also a gene-set average and should be similarly
# predictable. If it were NOT, the redundancy story would be specific
# to curated regulons rather than to gene-set averaging in general.
part_a <- list()
for (nm in names(tf_grid)) {
  tf_s <- scale_view(tf_grid[[nm]])
  r2 <- cv_r2_matrix(tf_s, X_pc, folds)
  r2 <- r2[is.finite(r2)]
  if (length(r2) == 0) next
  # Sensitivity: repeat with far fewer PCs. If redundancy is real it
  # should be visible at 50 PCs too. If R^2 only appears at 200, the
  # result is an artefact of predictor count and must not be reported
  # as redundancy. Nearly free -- the QR is the cost and it is small
  # at 50 predictors.
  r2_small <- cv_r2_matrix(tf_s, X_pc[, seq_len(min(PC_SENSITIVITY, ncol(X_pc))), drop = FALSE], folds)
  r2_small <- r2_small[is.finite(r2_small)]

  part_a[[nm]] <- tibble(
    view        = nm,
    is_permuted = grepl("edgeperm", nm),
    n_tf        = length(r2),
    median_cv_r2 = median(r2),
    q25_cv_r2    = quantile(r2, 0.25),
    q75_cv_r2    = quantile(r2, 0.75),
    frac_above_50 = mean(r2 > 0.50),
    frac_above_80 = mean(r2 > 0.80),
    median_cv_r2_50pc = if (length(r2_small)) median(r2_small) else NA_real_
  )
  cat(sprintf("  %-34s  n_TF=%3d  CV R2 = %.3f [%.3f, %.3f]  >0.8: %3.0f%%  (50 PCs: %.3f)\n",
              nm, length(r2), median(r2), quantile(r2, .25), quantile(r2, .75),
              100 * mean(r2 > 0.80),
              if (length(r2_small)) median(r2_small) else NA_real_))
  rm(tf_s, r2_small); gc(verbose = FALSE)
}
part_a_tbl <- bind_rows(part_a)
stopifnot(nrow(part_a_tbl) > 0)

# NEGATIVE CONTROL. Shuffle sample labels of one real view and repeat.
# If the CV estimator is honest, R^2 must collapse to ~0. Without this
# check, a high median R^2 could reflect a bug in the fold logic rather
# than genuine redundancy -- and the entire mechanistic claim would
# rest on that bug.
# Note on what to expect. Out-of-sample R^2 for unpredictable data is
# not 0 but slightly NEGATIVE, by roughly -(n_pc / n_train) -- a model
# fitted on noise predicts worse than the training mean. A value near
# -0.2 here is therefore the correct answer, not a problem. What would
# indicate a bug is a clearly POSITIVE value, which is what the check
# below tests for.
cat("\n  NEGATIVE CONTROL (sample labels shuffled; expect <= 0,\n")
cat("  mildly negative is normal; clearly positive would mean leakage):\n")
ctrl_view <- if ("dorothea_AC__viper" %in% names(tf_grid)) "dorothea_AC__viper" else names(tf_grid)[1]
set.seed(99001)
tf_shuf <- scale_view(tf_grid[[ctrl_view]])[sample(n_samples), , drop = FALSE]
r2_shuf <- cv_r2_matrix(tf_shuf, X_pc, folds)
r2_shuf <- r2_shuf[is.finite(r2_shuf)]
cat(sprintf("  %-34s  median CV R2 = %.4f  (expected ~0)\n",
            paste0(ctrl_view, " [shuffled]"), median(r2_shuf)))
cat(sprintf("  (for reference, -n_pc/n_train = %.3f)\n",
            -n_pc / (n_samples * (N_FOLDS - 1) / N_FOLDS)))
if (median(r2_shuf) > 0.10) {
  warning("Shuffled-control CV R2 is ", round(median(r2_shuf), 3),
          ", which should be ~0. The cross-validation is leaking. ",
          "DO NOT report Part A until this is resolved.")
}
rm(tf_shuf, pca); gc(verbose = FALSE)

part_a_tbl <- bind_rows(
  part_a_tbl,
  tibble(view = paste0(ctrl_view, "__SHUFFLED_CONTROL"), is_permuted = NA,
         n_tf = length(r2_shuf), median_cv_r2 = median(r2_shuf),
         q25_cv_r2 = quantile(r2_shuf, .25), q75_cv_r2 = quantile(r2_shuf, .75),
         frac_above_50 = mean(r2_shuf > .5), frac_above_80 = mean(r2_shuf > .8),
         median_cv_r2_50pc = NA_real_)
)
write_csv(part_a_tbl, here::here("results/tables/table_tf_expression_redundancy.csv"))

# ==================================================================
# PART B -- GEOMETRY-LEVEL REDUNDANCY
# ==================================================================
cat("\n=== PART B: how far does fusion move the geometry SNF acts on? ===\n\n")

D_expr <- dist2_cached(expr_scaled)
verify_distance_cache(expr_scaled, D_expr, label = "expression")
W_expr <- affinity_from_dist(D_expr, K = SNF_K, sigma = SNF_SIGMA)

lower <- lower.tri(D_expr)
d_expr_vec <- D_expr[lower]
w_expr_vec <- W_expr[lower]

# Expression-only reference partitions -- the counterfactual "what
# would we have got without the TF view at all".
cl_expr <- lapply(K_VALUES, function(k) SNFtool::spectralClustering(W_expr, K = k))
names(cl_expr) <- as.character(K_VALUES)

# ------------------------------------------------------------------
# PART B IS SPLIT IN TWO, because its two questions have wildly
# different costs and there is no reason to pay the high one 29 times.
#
#   B1. "How similar are the two views' geometries?" needs only the two
#       distance matrices. No SNF. About five seconds per view, so it
#       runs on ALL views and every number in it is complete.
#
#   B2. "How far does FUSION move the partition?" needs an actual SNF
#       run: 210 GFLOP on this machine's reference BLAS, which script
#       07's own log prices at roughly 15-20 minutes per view. It
#       therefore runs on a small pre-specified set only.
#
# This is not a shortcut. B1 answers the redundancy question for the
# whole grid; B2 confirms that the redundancy propagates through fusion,
# and confirming that on five views including two randomised controls is
# sufficient -- doing it on 29 would cost seven hours to restate a
# result already visible in the first five.
# ------------------------------------------------------------------

# ---------- B1: distance-level, ALL views, cheap ----------
cat("--- B1: geometry similarity, all views (no SNF; ~5 s/view) ---\n")
b1 <- list()
for (nm in names(tf_grid)) {
  tf_s <- scale_view(tf_grid[[nm]])
  D_tf <- dist2_cached(tf_s)
  b1[[nm]] <- tibble(
    view = nm, is_permuted = grepl("edgeperm", nm),
    rho_dist_expr_vs_tf = suppressWarnings(
      cor(d_expr_vec, D_tf[lower], method = "spearman")))
  cat(sprintf("  %-34s  rho(D_expr, D_TF) = %.3f\n",
              nm, b1[[nm]]$rho_dist_expr_vs_tf))
  rm(tf_s, D_tf); gc(verbose = FALSE)
}
b1_tbl <- bind_rows(b1)
write_csv(b1_tbl, here::here("results/tables/table_geometry_redundancy_all_views.csv"))
cat("  B1 complete and saved for all ", nrow(b1_tbl), " views.\n\n", sep = "")

# ---------- B2: fusion-level, priority views only, expensive ----------
# CHECKPOINT. Each view is saved as it completes and skipped on a
# re-run, so an interrupt costs nothing already earned.
B_CKPT <- here::here("results/objects/ckpt_18_partB.rds")
part_b <- if (file.exists(B_CKPT)) readRDS(B_CKPT) else list()

# A checkpoint written under different settings would silently mix two
# analyses. Stamp the configuration and discard the file if it changed.
b_stamp <- paste(SNF_K, SNF_SIGMA, SNF_T, paste(K_VALUES, collapse = "-"),
                 length(common_samples), sep = "|")
if (length(part_b) > 0 && !identical(attr(part_b, "stamp"), b_stamp)) {
  cat("  checkpoint was written under different settings -- discarding it\n")
  part_b <- list()
}
if (length(part_b) > 0)
  cat("  resuming: ", length(part_b), " view(s) already done\n", sep = "")

view_order <- intersect(GEOMETRY_PRIORITY, names(tf_grid))
if (RUN_ALL_VIEWS_IN_B2)
  view_order <- c(view_order, setdiff(names(tf_grid), view_order))

cat("--- B2: fusion effect, ", length(view_order), " view(s), ~15-20 min each ---\n",
    sep = "")
cat("    (SNF at n=", length(common_samples), ", t=", SNF_T,
    " is ~210 GFLOP; this is the slow part)\n", sep = "")

t_b0 <- Sys.time(); done_this_run <- 0L
for (i in seq_along(view_order)) {
  nm <- view_order[i]
  if (!is.null(part_b[[nm]])) next
  tf_s   <- scale_view(tf_grid[[nm]])
  D_tf   <- dist2_cached(tf_s)
  W_tf   <- affinity_from_dist(D_tf, K = SNF_K, sigma = SNF_SIGMA)
  W_fus  <- SNFtool::SNF(list(W_expr, W_tf), K = SNF_K, t = SNF_T)

  rho_dist  <- suppressWarnings(cor(d_expr_vec, D_tf[lower], method = "spearman"))
  rho_fused <- suppressWarnings(cor(w_expr_vec, W_fus[lower], method = "spearman"))

  row <- tibble(view = nm, is_permuted = grepl("edgeperm", nm),
                rho_dist_expr_vs_tf = rho_dist,
                rho_affinity_fused_vs_expr = rho_fused)

  for (k in K_VALUES) {
    cl_f <- SNFtool::spectralClustering(W_fus, K = k)
    row[[paste0("ari_fused_vs_expronly_k", k)]] <-
      ari_vs_reference(cl_f, cl_expr[[as.character(k)]], keep = NULL)
  }
  part_b[[nm]] <- row
  attr(part_b, "stamp") <- b_stamp
  saveRDS(part_b, B_CKPT)
  done_this_run <- done_this_run + 1L
  cat(sprintf("  [%d/%d] %-34s  rho(W_fused,W_expr)=%.3f  ARI(fused,expr-only) k2=%.3f k3=%.3f\n",
              i, length(view_order), nm, rho_fused,
              row$ari_fused_vs_expronly_k2, row$ari_fused_vs_expronly_k3))
  flush.console()

  if (done_this_run == 1L) {
    per <- as.numeric(difftime(Sys.time(), t_b0, units = "mins"))
    left <- sum(vapply(view_order, function(x) is.null(part_b[[x]]), logical(1)))
    cat(sprintf("  [timing] %.1f min/view MEASURED -> ~%.0f min (%.1f h) remaining\n",
                per, per * left, per * left / 60))
    flush.console()
  }
  rm(tf_s, D_tf, W_tf, W_fus); gc(verbose = FALSE)
}

part_b_tbl <- bind_rows(part_b[view_order[view_order %in% names(part_b)]])
stopifnot(nrow(part_b_tbl) > 0)
cat("\n  B2 complete for ", nrow(part_b_tbl), " view(s).\n\n", sep = "")
write_csv(part_b_tbl, here::here("results/tables/table_geometry_redundancy.csv"))

# ==================================================================
# PART C -- DOES SNF FUSION BEAT NAIVE CONCATENATION?
# ==================================================================
cat("\n=== PART C: fusion machinery vs simply gluing features together ===\n")
cat("    (reads script 07 output; no recomputation)\n\n")

grid_path <- here::here("results/tables/table_benchmark_grid_k2_k3.csv")
part_c_tbl <- NULL
if (file.exists(grid_path)) {
  g <- read_csv(grid_path, show_col_types = FALSE)
  base_ari <- g %>% filter(arm == "baseline") %>% select(k, baseline_ari = ari_vs_PAM50)
  part_c_tbl <- g %>%
    filter(arm %in% c("SNF_fused", "naive_concatenation")) %>%
    select(k, view, arm, ari_vs_PAM50) %>%
    pivot_wider(names_from = arm, values_from = ari_vs_PAM50) %>%
    { if (!all(c("SNF_fused", "naive_concatenation") %in% names(.)) ||
          !is.numeric(.$SNF_fused) || !is.numeric(.$naive_concatenation))
        stop("Part C: pivot produced non-numeric columns, which means the ",
             "benchmark grid has more than one row per (k, view, arm). ",
             "Inspect table_benchmark_grid_k2_k3.csv before trusting this.")
      . } %>%
    left_join(base_ari, by = "k") %>%
    mutate(fused_minus_concat   = SNF_fused - naive_concatenation,
           fused_minus_baseline = SNF_fused - baseline_ari,
           concat_minus_baseline = naive_concatenation - baseline_ari)
  write_csv(part_c_tbl, here::here("results/tables/table_fusion_vs_concatenation.csv"))

  for (k in K_VALUES) {
    s <- part_c_tbl %>% filter(k == !!k)
    if (nrow(s) == 0) next
    cat(sprintf("  k=%d  (n=%d views)\n", k, nrow(s)))
    cat(sprintf("    baseline ARI                     = %.4f\n", s$baseline_ari[1]))
    cat(sprintf("    SNF_fused    - concatenation     : median %+.4f  [%+.4f, %+.4f]\n",
                median(s$fused_minus_concat, na.rm = TRUE),
                min(s$fused_minus_concat, na.rm = TRUE),
                max(s$fused_minus_concat, na.rm = TRUE)))
    cat(sprintf("    SNF_fused    - baseline          : median %+.4f\n",
                median(s$fused_minus_baseline, na.rm = TRUE)))
    cat(sprintf("    concatenation- baseline          : median %+.4f\n",
                median(s$concat_minus_baseline, na.rm = TRUE)))
    cat(sprintf("    views where fusion beats concat  : %d / %d\n\n",
                sum(s$fused_minus_concat > 0, na.rm = TRUE), nrow(s)))
  }
} else {
  cat("  SKIPPED: ", grid_path, " not found.\n\n")
}

# ==================================================================
# SUMMARY
# ==================================================================
real_a <- part_a_tbl %>% filter(!is.na(is_permuted), !is_permuted)
perm_a <- part_a_tbl %>% filter(!is.na(is_permuted), is_permuted)
real_b <- part_b_tbl %>% filter(!is_permuted)

cat("\n============================================================\n")
cat("SUMMARY -- MECHANISM (descriptive, post hoc)\n")
cat("============================================================\n")
cat(sprintf("\nA. FEATURE LEVEL. Median cross-validated R^2 for predicting TF\n"))
cat(sprintf("   activity from %d expression PCs:\n", n_pc))
cat(sprintf("     real regulons          %.3f  (range %.3f - %.3f across views)\n",
            median(real_a$median_cv_r2), min(real_a$median_cv_r2), max(real_a$median_cv_r2)))
if (nrow(perm_a) > 0)
  cat(sprintf("     edge-permuted regulons %.3f\n", median(perm_a$median_cv_r2)))
cat(sprintf("     shuffled control       %.3f  (sanity check; must be <= 0)\n",
            median(r2_shuf)))
cat(sprintf("     TFs with R^2 > 0.8:    %.0f%% of real-regulon TFs\n",
            100 * median(real_a$frac_above_80)))

cat(sprintf("\nB. GEOMETRY LEVEL (real regulons):\n"))
cat(sprintf("     Spearman rho(D_expression, D_TF)        median %.3f  (ALL %d views, B1)\n",
            median(b1_tbl$rho_dist_expr_vs_tf[!b1_tbl$is_permuted], na.rm = TRUE),
            nrow(b1_tbl)))
cat(sprintf("     Spearman rho(W_fused, W_expression)     median %.3f\n",
            median(real_b$rho_affinity_fused_vs_expr, na.rm = TRUE)))
for (k in K_VALUES) {
  col <- paste0("ari_fused_vs_expronly_k", k)
  cat(sprintf("     ARI(fused partition, expression-only)   k=%d  median %.3f\n",
              k, median(real_b[[col]], na.rm = TRUE)))
}

cat("\nHOW TO READ THIS\n")
cat("----------------\n")
cat("High Part-A R^2 means the TF view is largely a linear re-encoding of\n")
cat("expression features the model already has. High Part-B ARI means the\n")
cat("fused partition is close to the one expression alone would have\n")
cat("produced. Together they convert the null of scripts 07 and 14 from an\n")
cat("observation into an explanation: there was little independent\n")
cat("information available to add, so no method operating on this feature\n")
cat("space could have added much.\n\n")
cat("If instead Part-A R^2 is LOW and Part-B ARI is HIGH, the conclusion\n")
cat("changes: independent information existed and SNF discarded it, which\n")
cat("is a criticism of the fusion method rather than of the TF view. Both\n")
cat("outcomes are reportable; they are reported differently.\n")

# ------------------------------------------------------------------
# PRODUCED OUTPUTS
# ------------------------------------------------------------------
cat("\n=== PRODUCED OUTPUTS ===\n")
outs <- c("results/tables/table_tf_expression_redundancy.csv",
          "results/tables/table_geometry_redundancy_all_views.csv",
          "results/tables/table_geometry_redundancy.csv",
          "results/tables/table_fusion_vs_concatenation.csv")
for (o in outs) {
  p <- here::here(o)
  cat(sprintf("  %-58s %s\n", o, if (file.exists(p)) "OK" else "MISSING"))
}
cat("\nScript 18 complete:", format(Sys.time()), "\n")
close_logger(log_con, script_name)
