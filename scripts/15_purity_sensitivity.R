# ============================================================
# 15_purity_sensitivity.R   (NEW)
# Is the cluster structure tumour biology, or tumour content?
# ============================================================
# WHY THIS SCRIPT EXISTS
# ------------------------------------------------------------
# Script 04 already reported, and then did nothing with, the single
# largest confound in this pipeline:
#
#     tumour purity vs PC1:  rho =  0.499   p = 2.8e-69
#     tumour purity vs PC2:  rho = -0.569   p = 4.8e-94
#
# Those are very large correlations against the two leading axes of the
# expression matrix that everything downstream is built on. A bulk
# RNA-seq sample is a mixture of tumour, stroma and immune cells, and
# the proportion varies enormously between samples. Any clustering of
# bulk expression can therefore separate samples by HOW MUCH TUMOUR IS
# IN THE TUBE rather than by what kind of tumour it is.
#
# A referee who reads the script-04 log and finds no follow-up will
# reasonably ask whether the "subtypes" are partly a stromal-content
# axis. That question must be answered in the manuscript, not left
# open, and the answer is valuable in EITHER direction:
#
#   clusters survive purity adjustment  -> a substantial credibility
#       gain, and a limitation pre-emptively closed
#   clusters do not survive             -> a finding in its own right,
#       and a considerably more interesting one than another
#       re-derivation of the ER+/ER- split
#
# The test is a purity-residualised re-clustering: regress purity out
# of every gene, re-select top-variance features on the residualised
# matrix, re-run the identical SNF procedure, and compare partitions by
# ARI. Feature RE-SELECTION matters and is easy to get wrong -- keeping
# the original top-2000 gene set would test only whether cluster
# ASSIGNMENTS move, while quietly assuming the FEATURE SELECTION itself
# was not purity-driven. Since variance is exactly what purity inflates,
# that assumption is the one most in need of testing, so both variants
# are run and reported separately.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SNFtool)
  library(mclust)
  library(survival)
})
source(here::here("R", "utils.R"))
source(here::here("R", "utils_benchmark.R"))

script_name <- "15_purity_sensitivity"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

SNF_K <- 20; SNF_SIGMA <- 0.5; SNF_T <- 20
K_VALUES <- c(2, 3)
N_TOP_GENES <- 2000

# ------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------
vst_all         <- readRDS(here::here("data/processed/vst_matrix_all_genes.rds"))
vst_top2000     <- readRDS(here::here("data/processed/vst_top2000_genes.rds"))
tf_activity     <- readRDS(here::here("data/processed/tf_activity_viper_AC_primary.rds"))
sample_metadata <- readRDS(here::here("data/processed/sample_metadata_matched.rds"))
final_labels    <- readRDS(here::here("results/objects/final_snf_cluster_labels.rds"))

purity_path <- here::here("results/tables/table_tumor_purity_merged.csv")
if (!file.exists(purity_path)) {
  stop("table_tumor_purity_merged.csv not found. Run the updated 04_quality_control.R ",
       "with data/external/tumor_purity_clean.csv in place first.")
}
purity_df <- read_csv(purity_path, show_col_types = FALSE) %>%
  dplyr::select(full_barcode, purity) %>%
  filter(!is.na(purity))

aligned <- align_samples(expr = vst_top2000, tf = tf_activity, meta = sample_metadata)
vst_top2000 <- aligned$expr; tf_activity <- aligned$tf; sample_metadata <- aligned$meta
samples_all <- colnames(vst_top2000)

purity_vec <- purity_df$purity[match(samples_all, purity_df$full_barcode)]
names(purity_vec) <- samples_all
has_purity <- !is.na(purity_vec)
cat("Samples with a purity estimate:", sum(has_purity), "of", length(samples_all), "\n")
if (sum(has_purity) < 300) stop("Too few samples with purity estimates for a meaningful test.")

# ------------------------------------------------------------------
# PART 1. Do the ORIGINAL clusters differ in purity?
#
# This is the descriptive question, and on its own it is NOT damning:
# real biological subtypes genuinely differ in stromal and immune
# content, so an association is expected. It becomes a problem only if
# the clustering is DRIVEN by purity, which Part 2 tests.
# ------------------------------------------------------------------
purity_assoc <- map_dfr(K_VALUES, function(k) {
  lab <- if (k == 2) final_labels$k2 else final_labels$k3
  cl <- lab$SNF_cluster[match(samples_all, lab$full_barcode)]
  ok <- has_purity & !is.na(cl)
  kw <- kruskal.test(purity_vec[ok] ~ factor(cl[ok]))
  eff <- map_dfr(sort(unique(cl[ok])), function(g) {
    tibble(cluster = g,
           median_purity = median(purity_vec[ok][cl[ok] == g]),
           n = sum(cl[ok] == g),
           cohens_d_vs_rest = cohens_d_one_vs_rest(purity_vec[ok], cl[ok], g))
  })
  bind_cols(tibble(k = k, kruskal_p = kw$p.value, kruskal_chisq = unname(kw$statistic)),
            eff %>% nest(cluster_detail = everything()))
})
purity_assoc_flat <- purity_assoc %>% unnest(cluster_detail)
cat("\n=== Purity by SNF cluster (original clustering) ===\n")
print(purity_assoc_flat, n = 20)
write_csv(purity_assoc_flat, here::here("results/tables/table_purity_by_cluster.csv"))

# ------------------------------------------------------------------
# PART 2. Purity-residualised re-clustering
#
# Per gene: expression ~ purity, keep residuals. This removes the
# linear purity component from every gene while preserving each gene's
# residual variation across samples. Linear adjustment is a deliberate,
# conservative choice: it is transparent and standard, but it will not
# remove non-linear or interaction-form purity effects, and that
# limitation is stated rather than glossed.
# ------------------------------------------------------------------
expr_sub <- vst_all[, samples_all[has_purity], drop = FALSE]
tf_sub   <- tf_activity[, samples_all[has_purity], drop = FALSE]
pur_sub  <- purity_vec[has_purity]

cat("\nResidualising", nrow(expr_sub), "genes on tumour purity ...\n")
design <- model.matrix(~ pur_sub)
qr_d <- qr(design)
resid_expr <- t(qr.resid(qr_d, t(expr_sub)))
dimnames(resid_expr) <- dimnames(expr_sub)

# TF activity is a function of expression, so it inherits the same
# confound and must be residualised too -- adjusting only the
# expression view would leave purity structure in the fused network by
# the back door.
resid_tf <- t(qr.resid(qr_d, t(tf_sub)))
dimnames(resid_tf) <- dimnames(tf_sub)

# Variant A: re-select top-variance genes ON THE RESIDUALISED matrix.
# This is the honest test, because it also asks whether the ORIGINAL
# feature selection was itself purity-driven.
resid_var <- apply(resid_expr, 1, var)
top_resid_genes <- names(sort(resid_var, decreasing = TRUE))[1:min(N_TOP_GENES, length(resid_var))]
expr_A <- resid_expr[top_resid_genes, , drop = FALSE]

# Variant B: keep the ORIGINAL top-2000 gene set, residualised. Isolates
# the effect of adjustment alone, holding features fixed.
orig_genes <- intersect(rownames(vst_top2000), rownames(resid_expr))
expr_B <- resid_expr[orig_genes, , drop = FALSE]

overlap_genes <- length(intersect(top_resid_genes, orig_genes))
cat("Feature-set overlap after residualisation:", overlap_genes, "of",
    length(top_resid_genes), sprintf(" (%.1f%%)\n", 100 * overlap_genes / length(top_resid_genes)))
cat("  -> a LOW overlap would itself indicate that the original variance-based\n",
    "     feature selection was substantially purity-driven.\n", sep = "")

scale_view <- function(mat) { s <- scale(t(mat)); s[is.na(s)] <- 0; s }

run_snf <- function(expr_mat, tf_mat, k) {
  W_e <- affinity_from_dist(dist2_cached(scale_view(expr_mat)), K = SNF_K, sigma = SNF_SIGMA)
  W_t <- affinity_from_dist(dist2_cached(scale_view(tf_mat)),   K = SNF_K, sigma = SNF_SIGMA)
  SNFtool::spectralClustering(SNFtool::SNF(list(W_e, W_t), K = SNF_K, t = SNF_T), K = k)
}

sensitivity_rows <- map_dfr(K_VALUES, function(k) {
  lab <- if (k == 2) final_labels$k2 else final_labels$k3
  orig_cl <- lab$cluster_int[match(colnames(expr_sub), lab$full_barcode)]

  cl_A <- run_snf(expr_A, resid_tf, k)
  cl_B <- run_snf(expr_B, resid_tf, k)

  # Residual purity association after adjustment: if this is still
  # significant, the linear adjustment did not fully remove the effect.
  kw_A <- kruskal.test(pur_sub ~ factor(cl_A))
  kw_B <- kruskal.test(pur_sub ~ factor(cl_B))

  tibble(
    k = k,
    variant = c("A_reselected_features", "B_original_features"),
    ari_vs_original = c(mclust::adjustedRandIndex(cl_A, orig_cl),
                        mclust::adjustedRandIndex(cl_B, orig_cl)),
    residual_purity_kruskal_p = c(kw_A$p.value, kw_B$p.value),
    n_samples = length(pur_sub)
  )
}) %>%
  mutate(
    interpretation = case_when(
      ari_vs_original >= 0.80 ~ "ROBUST: partition largely reproduced after purity adjustment",
      ari_vs_original >= 0.50 ~ "PARTIAL: partition materially shifts; purity is a substantial contributor",
      TRUE                     ~ "FRAGILE: partition largely dissolves; the original clusters are substantially a purity/stromal-content axis"
    )
  )

cat("\n=== Purity-residualised re-clustering ===\n")
print(sensitivity_rows, width = Inf)
write_csv(sensitivity_rows, here::here("results/tables/table_purity_sensitivity.csv"))

# ------------------------------------------------------------------
# PART 3. Purity as a competing explanation in the survival model
#
# If purity alone predicts survival as well as the cluster labels do,
# the clusters are not carrying independent prognostic information --
# a distinct question from whether the partition is stable, and one
# that a purely clustering-side sensitivity analysis cannot answer.
# ------------------------------------------------------------------
sm <- sample_metadata
sm$days_to_death         <- suppressWarnings(as.numeric(sm$days_to_death))
sm$days_to_last_followup <- suppressWarnings(as.numeric(sm$days_to_last_followup))
surv_years <- ifelse(tolower(sm$vital_status) == "dead",
                     sm$days_to_death, sm$days_to_last_followup) / 365.25
surv_status <- ifelse(tolower(sm$vital_status) == "dead", 1, 0)

surv_purity <- map_dfr(K_VALUES, function(k) {
  lab <- if (k == 2) final_labels$k2 else final_labels$k3
  cl <- lab$SNF_cluster[match(samples_all, lab$full_barcode)]
  d <- tibble(time = surv_years, status = surv_status, cl = factor(cl),
              purity = purity_vec, age = sm$age_years) %>%
    filter(is.finite(time), time > 0, !is.na(status), !is.na(cl), !is.na(purity), !is.na(age))
  if (nrow(d) < 50) return(tibble())
  f_cl   <- coxph(Surv(time, status) ~ cl + age, data = d)
  f_pur  <- coxph(Surv(time, status) ~ purity + age, data = d)
  f_both <- coxph(Surv(time, status) ~ cl + purity + age, data = d)
  tibble(
    k = k, n = nrow(d), events = sum(d$status),
    cindex_cluster_age = unname(summary(f_cl)$concordance[1]),
    cindex_purity_age  = unname(summary(f_pur)$concordance[1]),
    cindex_both        = unname(summary(f_both)$concordance[1]),
    lrt_cluster_over_purity_age = 2 * (as.numeric(logLik(f_both)) - as.numeric(logLik(f_pur))),
    lrt_cluster_over_purity_age_p = pchisq(
      2 * (as.numeric(logLik(f_both)) - as.numeric(logLik(f_pur))),
      df = length(unique(d$cl)) - 1, lower.tail = FALSE)
  )
})
cat("\n=== Cluster vs purity as prognostic factors (TCGA) ===\n")
print(surv_purity, width = Inf)
write_csv(surv_purity, here::here("results/tables/table_purity_vs_cluster_prognostic.csv"))

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
plot_df <- map_dfr(K_VALUES, function(k) {
  lab <- if (k == 2) final_labels$k2 else final_labels$k3
  tibble(k = paste0("k = ", k),
         cluster = lab$SNF_cluster[match(samples_all, lab$full_barcode)],
         purity = purity_vec)
}) %>% filter(!is.na(purity), !is.na(cluster))

p <- ggplot(plot_df, aes(x = cluster, y = purity, fill = cluster)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  facet_wrap(~ k, scales = "free_x") +
  labs(x = NULL, y = "Tumour purity",
       title = "Tumour purity by SNF cluster",
       subtitle = "A difference here is expected for real subtypes; see table_purity_sensitivity.csv for whether the partition is DRIVEN by purity") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggsave(here::here("results/figures/figure_purity_by_cluster.png"), p,
       width = 9, height = 5, dpi = 150)

writeLines(c(
  "TUMOUR-PURITY CONFOUND: what was tested and how to report it.",
  "======================================================================",
  "Script 04 reported very large correlations between ABSOLUTE tumour",
  "purity and the two leading expression principal components",
  "(rho = 0.50 with PC1; rho = -0.57 with PC2) and did not follow them",
  "up. Bulk RNA-seq measures a mixture of tumour, stroma and immune",
  "cells, so any clustering of bulk expression can separate samples by",
  "tumour CONTENT rather than tumour BIOLOGY. This script closes that",
  "gap.",
  "",
  "THREE SEPARATE QUESTIONS, THREE SEPARATE TABLES.",
  "1. table_purity_by_cluster.csv -- do clusters differ in purity?",
  "   An association alone is NOT evidence of artefact: real subtypes",
  "   genuinely differ in stromal and immune content.",
  "2. table_purity_sensitivity.csv -- does the partition SURVIVE having",
  "   purity regressed out of every gene? This is the decisive test.",
  "   Variant A also re-selects top-variance features on the",
  "   residualised matrix, and so additionally tests whether the",
  "   original variance-based feature selection was itself",
  "   purity-driven -- which the fixed-feature Variant B cannot detect.",
  "3. table_purity_vs_cluster_prognostic.csv -- do the clusters carry",
  "   prognostic information beyond purity alone?",
  "",
  "LIMITATION. Adjustment is LINEAR in purity. Non-linear purity",
  "effects, and interactions between purity and biology, are not",
  "removed. A cell-type deconvolution approach would address more of",
  "this but introduces its own reference-profile assumptions.",
  "",
  "HOW TO REPORT. Whichever way this comes out, report it in the main",
  "text rather than the supplement. If the clusters are robust, this is",
  "a strong pre-emptive answer to an obvious referee question. If they",
  "are not, that is a more interesting finding than another",
  "re-derivation of the ER+/ER- split, and should be framed as such."
), here::here("results/tables/NOTE_purity_sensitivity.txt"))

log_session_info(script_name, key_packages = c("SNFtool", "mclust", "survival"))
cat("\n✓ 15_purity_sensitivity.R complete.\n")

close_logger(log_con, script_name)
