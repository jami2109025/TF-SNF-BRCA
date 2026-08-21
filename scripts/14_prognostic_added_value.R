# ============================================================
# 14_prognostic_added_value.R   (NEW)
# Added value of TF-activity integration on a NON-CIRCULAR,
# OUTCOME-BASED endpoint in an independent cohort (METABRIC)
# ============================================================
# WHY THIS SCRIPT EXISTS
# ------------------------------------------------------------
# Script 07 measures added value as agreement with PAM50. That endpoint
# has a structural defect that no amount of extra bootstrapping can fix:
#
#     PAM50 labels are themselves derived from gene expression.
#     An expression-only clustering is therefore near-optimal for that
#     target BY CONSTRUCTION. The baseline arm is scored against a
#     label computed from its own input; the TF arm is not.
#
# A benchmark in which one arm plays at home is not a fair test, and a
# null result from it is close to uninformative. Reviewers will make
# exactly this point, and they will be right.
#
# This script measures added value on an endpoint that is not
# downstream of the expression matrix at all -- OVERALL SURVIVAL -- and
# does so in a COHORT THAT PLAYED NO PART IN DERIVING THE CLUSTERS.
# Neither arm has any structural advantage. Three defects of the
# original design are removed simultaneously:
#
#   circularity  -> the outcome is not computed from expression
#   power        -> METABRIC has ~1144 events vs ~150 in TCGA, and
#                   follow-up out to ~29 years vs a few years
#   relevance    -> discrimination of survival is the quantity a
#                   clinical reader actually cares about, whereas
#                   agreement with another expression-based classifier
#                   is of methodological interest only
#
# This is the analysis most likely to turn the thesis's central
# question from a null into a real finding in EITHER direction. If TF
# integration improves prognostic discrimination, that is a positive
# result the PAM50 endpoint was incapable of detecting. If it does not,
# the negative result is now supported on the endpoint that matters,
# and is far harder to dismiss.
#
# BOOTSTRAP DESIGN NOTE
# The TCGA-derived centroids are treated as FIXED, pre-trained objects
# and the bootstrap resamples METABRIC patients only. This is the
# correct frame for "how well does this fixed classifier discriminate
# in new data", and it is the frame a clinical reader assumes. It does
# NOT propagate uncertainty in the discovery-cohort clustering itself;
# a nested bootstrap that re-derived clusters per replicate would, at
# very large computational cost. Stated as a limitation, not hidden.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(org.Hs.eg.db)
  library(survival)
  library(SNFtool)
})
source(here::here("R", "utils.R"))
source(here::here("R", "utils_benchmark.R"))

script_name <- "14_prognostic_added_value"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

# ------------------------------------------------------------------
# Self-contained parallel helpers (deliberately NOT in utils_benchmark.R)
#
# Script 07 is running concurrently and its workers re-source
# R/utils_benchmark.R when each new chunk starts. Editing that shared
# file now would inject untested code into a job with ~30 hours left to
# run, so these helpers are duplicated here instead. Duplication is the
# cheaper risk.
#
# What they protect against, both observed in script 07 today:
#   - makeCluster() blocking forever on a localhost socket (no error,
#     0% CPU, indistinguishable from "still working")
#   - a worker-side failure surfacing only after the full job has been
#     dispatched, ~15 minutes in
# ------------------------------------------------------------------
make_cluster_checked_14 <- function(n, timeout_sec = 600) {
  cl <- tryCatch({
    setTimeLimit(elapsed = timeout_sec, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
    parallel::makeCluster(n)
  }, error = function(e) {
    cat("  [parallel] makeCluster failed/timed out:", conditionMessage(e), "\n"); NULL
  })
  if (is.null(cl)) return(NULL)
  ok <- tryCatch({
    r <- parallel::clusterEvalQ(cl, 1L + 1L)
    length(r) == n && all(unlist(r) == 2L)
  }, error = function(e) FALSE)
  if (!ok) {
    cat("  [parallel] round-trip test FAILED\n")
    try(parallel::stopCluster(cl), silent = TRUE); return(NULL)
  }
  cat("  [parallel] cluster of", n, "workers verified.\n"); cl
}

par_or_serial_14 <- function(cl, X, FUN, label = "") {
  if (!is.null(cl)) return(parallel::parLapply(cl, X, FUN))
  cat("  [serial] running", length(X), "task(s) sequentially with progress.\n")
  out <- vector("list", length(X))
  for (i in seq_along(X)) {
    t0 <- Sys.time(); out[[i]] <- FUN(X[[i]])
    cat(sprintf("    [%s] %d/%d (%.1f min)\n", label, i, length(X),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    utils::flush.console()
  }
  out
}

# ------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------
# Replicates for the METABRIC bootstrap.
#
# Set to 2000 when the grid was expected to hold ~6 views. Script 05
# actually produced 29, and this bootstrap is SERIAL and scales as
# N_BOOT x n_views: at 2000 it is ~2 hours of Cox refitting on top of
# the ~1 hour of SNF above. 1000 halves that at negligible statistical
# cost -- the percentile CIs and the TOST read the 2.5/97.5 and 5/95
# quantiles, and 1000 replicates put 25 and 50 observations in those
# tails respectively, which is ample. (The floor to avoid is ~200, where
# the tails become unstable; 1000 is nowhere near it.)
N_BOOT_METABRIC <- 1000
K_VALUES        <- c(2, 3)
SNF_K <- 20; SNF_SIGMA <- 0.5; SNF_T <- 20
MIN_COVERAGE_PCT <- 70
SESOI_CINDEX <- 0.02   # pre-specified; see justification in the note at end

# ------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------
metabric_expr_path <- here::here("data/external/metabric_expression.txt")
metabric_clin_path <- here::here("data/external/metabric_clinical.txt")
if (!file.exists(metabric_expr_path) || !file.exists(metabric_clin_path)) {
  stop("METABRIC files not found under data/external/. See script 13 header.")
}

vst_top2000          <- readRDS(here::here("data/processed/vst_top2000_genes.rds"))
vst_top5000          <- readRDS(here::here("data/processed/vst_top5000_genes.rds"))
vst_matrix_all_genes <- readRDS(here::here("data/processed/vst_matrix_all_genes.rds"))
tf_grid              <- readRDS(here::here("data/processed/tf_activity_grid.rds"))
sample_metadata      <- readRDS(here::here("data/processed/sample_metadata_matched.rds"))

common_samples <- Reduce(intersect, c(
  list(colnames(vst_top2000), colnames(vst_top5000), rownames(sample_metadata)),
  lapply(tf_grid, colnames)
))
vst_top2000 <- vst_top2000[, common_samples, drop = FALSE]
vst_top5000 <- vst_top5000[, common_samples, drop = FALSE]
tf_grid     <- lapply(tf_grid, function(m) m[, common_samples, drop = FALSE])
cat("TCGA discovery samples used:", length(common_samples), "\n")

# ------------------------------------------------------------------
# STEP 1. Derive TCGA cluster labels under every arm
# ------------------------------------------------------------------
scale_view <- function(mat) { s <- scale(t(mat)); s[is.na(s)] <- 0; s }
expr_scaled <- scale_view(vst_top2000)
D_expr <- dist2_cached(expr_scaled)
verify_distance_cache(expr_scaled, D_expr, label = "expression")

# LOOP ORDER AND PARALLELISM. An earlier version of this block nested
# the view loop inside a k loop, which rebuilt every fused network once
# per k value even though SNF does not depend on k -- only the final
# spectralClustering(K = k) does. With a full grid and SNF measured at
# ~700 s per call on this cohort, that mistake alone costs several
# hours. Networks are now built once per view, clustered at every k,
# and the whole view loop is parallelised.
arm_labels <- list()

W_expr <- affinity_from_dist(D_expr, K = SNF_K, sigma = SNF_SIGMA)
for (k in K_VALUES) {
  arm_labels[[paste0("expression_only|k", k)]] <-
    setNames(SNFtool::spectralClustering(W_expr, K = k), common_samples)
}

# CORE BUDGET. Script 07 is running with 7 workers on an 8-core machine
# and is on the critical path. Requesting 7 more here would put 14
# compute processes on 8 cores: both jobs would thrash, and 07 -- the
# one that cannot be restarted in the time remaining -- would suffer
# most. Three workers leaves 07 the bulk of the machine while still
# cutting this stage from ~6 h serial to ~2 h. Raise it only once 07
# has finished.
# Script 07 has finished, so the machine is free -- but it still has
# only 8 GB. With payload passing each worker costs ~310 MB and the
# master holds the METABRIC expression matrix (~320 MB) plus the VST
# matrices. Five workers is ~1.6 GB of workers against ~3.4 GB free,
# which leaves real headroom; seven would leave essentially none.
N_CORES_14 <- 5
n_cores_14 <- max(1, min(N_CORES_14, length(tf_grid)))
cat("Deriving TCGA labels for", length(tf_grid), "views on", n_cores_14, "cores ...\n")

# MEMORY. tf_grid is 73 MB. Exporting it whole put a full copy in every
# worker -- 511 MB across seven, on a machine with 3.4 GB free. Each
# view's matrix is now passed as a parLapply ARGUMENT, so a worker only
# ever holds the views it is actually processing (~3 MB each).
view_label_fun <- function(payload) {
  nm <- payload$nm
  s <- scale(t(payload$mat)); s[is.na(s)] <- 0
  D_tf <- SNFtool::dist2(as.matrix(s), as.matrix(s))
  W_tf <- affinity_from_dist(D_tf, K = SNF_K, sigma = SNF_SIGMA)
  W_fused <- SNFtool::SNF(list(W_expr, W_tf), K = SNF_K, t = SNF_T)
  out <- lapply(K_VALUES, function(k) SNFtool::spectralClustering(W_fused, K = k))
  names(out) <- paste0("k", K_VALUES)
  out
}

cl14 <- make_cluster_checked_14(n_cores_14)
if (!is.null(cl14)) {
parallel::clusterEvalQ(cl14, {
  suppressPackageStartupMessages(library(SNFtool))
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
  }
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
})
parallel::clusterExport(
  cl14,
  varlist = c("W_expr", "SNF_K", "SNF_SIGMA", "SNF_T", "K_VALUES",
              "AFFINITY_EPS", "affinity_from_dist"),
  envir = environment()
)
}

# SMOKE TEST -- one view before committing to all 29.
cat("  [smoke] deriving labels for one view first ...\n"); utils::flush.console()
view_payloads_14 <- lapply(names(tf_grid), function(nm) list(nm = nm, mat = tf_grid[[nm]]))
names(view_payloads_14) <- names(tf_grid)

smoke14 <- tryCatch(par_or_serial_14(cl14, view_payloads_14[1], view_label_fun,
                                     label = "smoke")[[1]],
                    error = function(e) e)
if (inherits(smoke14, "error")) {
  if (!is.null(cl14)) try(parallel::stopCluster(cl14), silent = TRUE)
  stop("SMOKE TEST FAILED before deriving TCGA labels.\n  Message: ",
       conditionMessage(smoke14))
}
cat("  [smoke] PASS -- cluster sizes k=2:",
    paste(table(smoke14$k2), collapse = "/"), " k=3:",
    paste(table(smoke14$k3), collapse = "/"), "\n")
utils::flush.console()

# RETRY. parLapply signals if any worker dies, and a death 45 minutes
# into a 50-minute stage would otherwise discard the whole stage. Two
# attempts on a fresh cluster, then serial.
view_labels <- NULL
for (attempt in 1:2) {
  view_labels <- tryCatch(
    par_or_serial_14(cl14, view_payloads_14, view_label_fun, label = "labels"),
    error = function(e) e)
  if (!inherits(view_labels, "error")) break
  cat(sprintf("  [retry %d/2] label derivation failed: %s\n",
              attempt, conditionMessage(view_labels)))
  utils::flush.console()
  try(parallel::stopCluster(cl14), silent = TRUE)
  gc(verbose = FALSE, full = TRUE); Sys.sleep(15)
  cl14 <- make_cluster_checked_14(n_cores_14)
  if (!is.null(cl14)) {
    parallel::clusterEvalQ(cl14, {
      suppressPackageStartupMessages(library(SNFtool))
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
      }
      Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
    })
    parallel::clusterExport(cl14,
      varlist = c("W_expr", "SNF_K", "SNF_SIGMA", "SNF_T", "K_VALUES",
                  "AFFINITY_EPS", "affinity_from_dist"),
      envir = environment())
  }
}
if (inherits(view_labels, "error")) {
  cat("  [fallback] running label derivation serially.\n")
  view_labels <- lapply(view_payloads_14, view_label_fun)
}
if (!is.null(cl14)) try(parallel::stopCluster(cl14), silent = TRUE)
names(view_labels) <- names(tf_grid)

for (nm in names(view_labels)) {
  for (k in K_VALUES) {
    arm_labels[[paste0(nm, "|k", k)]] <-
      setNames(view_labels[[nm]][[paste0("k", k)]], common_samples)
  }
}
cat("TCGA labels derived across", length(arm_labels), "arm-by-k combinations\n")

# ------------------------------------------------------------------
# STEP 2. TCGA centroids in the top-5000 gene-symbol space
# (identical projection machinery to script 13, so the two external
# validations remain directly comparable)
# ------------------------------------------------------------------
mapped_5000 <- map_to_symbol_dedup(vst_top5000, variance_reference = vst_matrix_all_genes,
                                   org_db = org.Hs.eg.db)
tcga_z <- t(scale(t(mapped_5000$matrix)))
# vst_matrix_all_genes (~200 MB) was needed only as the variance
# reference for the symbol mapping above. On an 8 GB machine it must not
# be held for the rest of the script.
rm(vst_matrix_all_genes, mapped_5000); gc(verbose = FALSE, full = TRUE)
cat("TCGA centroid feature space:", nrow(tcga_z), "gene symbols\n")

build_centroids <- function(labels_vec, expr_z) {
  ids <- names(labels_vec)
  expr_z <- expr_z[, ids, drop = FALSE]
  sapply(split(ids, paste0("SNF_C", labels_vec)), function(s) {
    rowMeans(expr_z[, s, drop = FALSE])
  })
}

# ------------------------------------------------------------------
# STEP 3. Load and z-score METABRIC (gene-level variance filtering
# BEFORE scaling, as fixed in script 13)
# ------------------------------------------------------------------
metabric_raw <- read.delim(metabric_expr_path, check.names = FALSE)
gene_col <- intersect(c("Hugo_Symbol", "GENE_SYMBOL", "gene_symbol"), colnames(metabric_raw))[1]
if (is.na(gene_col)) stop("No gene-symbol column found in METABRIC expression file.")
known_id_cols <- c("Entrez_Gene_Id", "ENTREZ_GENE_ID", gene_col)

metabric_expr <- metabric_raw %>%
  dplyr::select(-any_of(setdiff(known_id_cols, gene_col))) %>%
  dplyr::rename(SYMBOL = all_of(gene_col)) %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()
rm(metabric_raw); gc(verbose = FALSE, full = TRUE)
cat("METABRIC expression matrix:", dim(metabric_expr), "\n")

common_genes <- intersect(rownames(tcga_z), rownames(metabric_expr))
coverage_pct <- 100 * length(common_genes) / nrow(tcga_z)
cat("Gene coverage:", length(common_genes), "/", nrow(tcga_z),
    sprintf(" (%.1f%%)\n", coverage_pct))
if (coverage_pct < MIN_COVERAGE_PCT) stop("Gene coverage below ", MIN_COVERAGE_PCT, "%.")

metabric_common <- metabric_expr[common_genes, ]
gene_var <- apply(metabric_common, 1, var, na.rm = TRUE)
drop_genes <- names(gene_var)[is.na(gene_var) | gene_var == 0]
if (length(drop_genes)) {
  metabric_common <- metabric_common[setdiff(rownames(metabric_common), drop_genes), ]
  common_genes <- setdiff(common_genes, drop_genes)
}
metabric_z <- t(scale(t(metabric_common)))
bad_post <- rownames(metabric_z)[rowSums(is.na(metabric_z)) > 0]
if (length(bad_post)) {
  metabric_z <- metabric_z[setdiff(rownames(metabric_z), bad_post), ]
  common_genes <- setdiff(common_genes, bad_post)
}
cat("Final common gene set:", length(common_genes), "\n")

# ------------------------------------------------------------------
# STEP 4. METABRIC clinical / survival
# ------------------------------------------------------------------
metabric_clin <- read.delim(metabric_clin_path, check.names = FALSE, comment.char = "#")

id_col <- intersect(c("PATIENT_ID", "SAMPLE_ID", "Patient.ID"), colnames(metabric_clin))[1]
if (is.na(id_col)) stop("No patient/sample ID column found in METABRIC clinical file.")

os_months_col <- intersect(c("OS_MONTHS", "os_months"), colnames(metabric_clin))[1]
os_status_col <- intersect(c("OS_STATUS", "os_status"), colnames(metabric_clin))[1]
if (is.na(os_months_col) || is.na(os_status_col)) {
  stop("OS_MONTHS / OS_STATUS not found in METABRIC clinical file.")
}

age_col <- pick_best_covariate(c("AGE_AT_DIAGNOSIS", "Age.at.Diagnosis"), metabric_clin)
er_col  <- pick_best_covariate(c("ER_IHC", "ER_STATUS", "ER_status"), metabric_clin)
if (is.na(age_col) || is.na(er_col)) {
  stop("Could not locate an age and/or ER column in the METABRIC clinical file. ",
       "Adjusted analysis is the whole point of this script -- an unadjusted ",
       "comparison would reproduce the ER-driven result of script 13 and must ",
       "not be substituted silently. Columns available: ",
       paste(colnames(metabric_clin), collapse = ", "))
}
cat("METABRIC covariates selected by completeness -- age:", age_col, "| ER:", er_col, "\n")

clin <- metabric_clin %>%
  transmute(
    sample_id = .data[[id_col]],
    os_months = suppressWarnings(as.numeric(.data[[os_months_col]])),
    os_status_raw = as.character(.data[[os_status_col]]),
    age_covariate = suppressWarnings(as.numeric(.data[[age_col]])),
    er_covariate = as.character(.data[[er_col]])
  ) %>%
  mutate(
    os_status = case_when(
      grepl("^1|DECEASED|DEAD", os_status_raw, ignore.case = TRUE) ~ 1,
      grepl("^0|LIVING|ALIVE", os_status_raw, ignore.case = TRUE) ~ 0,
      TRUE ~ NA_real_
    ),
    os_years = os_months / 12,
    # cBioPortal's METABRIC export contains the misspelling "Positve";
    # left unhandled it silently creates a third ER level and changes
    # the model's degrees of freedom.
    er_covariate = case_when(
      grepl("^posit", er_covariate, ignore.case = TRUE) ~ "Positive",
      grepl("^negat", er_covariate, ignore.case = TRUE) ~ "Negative",
      TRUE ~ NA_character_
    )
  )

cat("METABRIC OS status resolved:\n"); print(table(clin$os_status, useNA = "ifany"))

# --- EXPORT FOR SCRIPT 19 -------------------------------------------
# scripts/19_nested_bootstrap_sensitivity.R needs this exact clinical
# frame. Re-deriving the recoding there would risk the two analyses
# silently diverging on OS-status or ER coding -- note the cBioPortal
# misspelling "Positve" handled above -- and a divergence there is
# precisely the confound a sensitivity analysis must not introduce.
saveRDS(clin, here::here("results/objects/metabric_clin_prepared.rds"))
cat("Exported prepared clinical frame for script 19:",
    nrow(clin), "rows ->", "results/objects/metabric_clin_prepared.rds\n")
# --------------------------------------------------------------------

# ------------------------------------------------------------------
# STEP 5. Project METABRIC under every arm, score prognostic
#         discrimination, and bootstrap the PAIRED difference
# ------------------------------------------------------------------
classify_to_centroids <- function(mz, centroids, genes) {
  cmat <- centroids[genes, , drop = FALSE]
  cor_matrix <- cor(mz[genes, ], cmat, method = "pearson")
  best_idx <- max.col(cor_matrix, ties.method = "first")
  setNames(colnames(cmat)[best_idx], rownames(cor_matrix))
}

assignments <- list()
for (arm in names(arm_labels)) {
  cents <- build_centroids(arm_labels[[arm]], tcga_z)
  assignments[[arm]] <- classify_to_centroids(metabric_z, cents, common_genes)
}
cat("METABRIC projected under", length(assignments), "arms\n")

# Analysis frame: one row per METABRIC patient, complete on outcome
# and both covariates, so every arm is scored on the IDENTICAL patient
# set. Letting the analysable set vary by arm would confound the
# comparison with differential missingness.
base_df <- tibble(sample_id = colnames(metabric_z)) %>%
  left_join(clin, by = "sample_id") %>%
  filter(is.finite(os_years), os_years > 0, !is.na(os_status),
         !is.na(age_covariate), !is.na(er_covariate))

cat("METABRIC analysis set (complete on outcome + age + ER):", nrow(base_df),
    "patients,", sum(base_df$os_status == 1), "events\n")
cat("Median follow-up (years):", round(median(base_df$os_years), 2), "\n")

score_arm <- function(cluster_vec, df, rows) {
  d <- df[rows, , drop = FALSE]
  d$cl <- factor(cluster_vec[d$sample_id])
  if (length(unique(d$cl)) < 2) return(c(cindex_adj = NA_real_, cindex_cl = NA_real_, lrt = NA_real_))
  fit_cov  <- try(coxph(Surv(os_years, os_status) ~ age_covariate + er_covariate, data = d), silent = TRUE)
  fit_full <- try(coxph(Surv(os_years, os_status) ~ cl + age_covariate + er_covariate, data = d), silent = TRUE)
  fit_cl   <- try(coxph(Surv(os_years, os_status) ~ cl, data = d), silent = TRUE)
  if (inherits(fit_full, "try-error") || inherits(fit_cov, "try-error")) {
    return(c(cindex_adj = NA_real_, cindex_cl = NA_real_, lrt = NA_real_))
  }
  c(
    cindex_adj = unname(summary(fit_full)$concordance[1]),
    cindex_cl  = if (inherits(fit_cl, "try-error")) NA_real_ else unname(summary(fit_cl)$concordance[1]),
    # LRT of the cluster term OVER age + ER: the incremental prognostic
    # information the clustering carries beyond routine clinical
    # variables. A cluster solution that merely re-encodes ER status
    # scores ~0 here, which is exactly the discrimination the TCGA-only
    # analysis could not make.
    lrt = 2 * (as.numeric(logLik(fit_full)) - as.numeric(logLik(fit_cov)))
  )
}

all_rows <- seq_len(nrow(base_df))
point_scores <- map_dfr(names(assignments), function(arm) {
  s <- score_arm(assignments[[arm]], base_df, all_rows)
  parts <- strsplit(arm, "\\|")[[1]]
  tibble(view = parts[1], k = as.integer(sub("^k", "", parts[2])),
         cindex_adjusted = s[["cindex_adj"]], cindex_cluster_only = s[["cindex_cl"]],
         lrt_over_age_er = s[["lrt"]])
})
write_csv(point_scores, here::here("results/tables/table_metabric_prognostic_point_scores.csv"))

cat("\n=== METABRIC prognostic discrimination by arm (point estimates) ===\n")
print(point_scores %>% arrange(k, desc(cindex_adjusted)), n = 40)

# Paired bootstrap over METABRIC patients. All arms are evaluated on
# the SAME resampled patients within a replicate, so the difference is
# a paired contrast.
cat("\nRunning paired METABRIC bootstrap (n =", N_BOOT_METABRIC, ") ...\n")
t0 <- Sys.time()

# BATCHED WITH RESUME. This loop is serial and takes about an hour; a
# crash at minute 55 would otherwise discard all of it. Batches of 100
# are written to disk as they complete and skipped on re-run.
MB_BATCH <- 100
mb_starts <- seq(1, N_BOOT_METABRIC, by = MB_BATCH)
mb_files  <- character(0)

for (mi in seq_along(mb_starts)) {
  lo <- mb_starts[mi]; hi <- min(lo + MB_BATCH - 1, N_BOOT_METABRIC)
  mbf <- here::here("results/objects",
                    sprintf("metabric_boot_%04d_%04d.rds", lo, hi))
  mb_files <- c(mb_files, mbf)
  if (file.exists(mbf)) { cat(sprintf("  batch %d-%d already on disk, skipping\n", lo, hi)); next }

  t_mb <- Sys.time()
  boot_rows <- vector("list", hi - lo + 1)
for (b in lo:hi) {
  set.seed(PIPELINE_SEED + 6000 + b)
  rows <- sample(all_rows, length(all_rows), replace = TRUE)
  rep_out <- list()
  for (k in K_VALUES) {
    base_arm <- paste0("expression_only|k", k)
    base_s <- score_arm(assignments[[base_arm]], base_df, rows)
    for (nm in names(tf_grid)) {
      s <- score_arm(assignments[[paste0(nm, "|k", k)]], base_df, rows)
      rep_out[[length(rep_out) + 1]] <- data.frame(
        b = b, k = k, view = nm,
        d_cindex_adj = s[["cindex_adj"]] - base_s[["cindex_adj"]],
        d_cindex_cl  = s[["cindex_cl"]]  - base_s[["cindex_cl"]],
        d_lrt        = s[["lrt"]]        - base_s[["lrt"]]
      )
    }
  }
  boot_rows[[b - lo + 1]] <- do.call(rbind, rep_out)
}
  saveRDS(dplyr::bind_rows(boot_rows), mbf)
  cat(sprintf("  batch %d-%d done (%.1f min) | %d/%d replicates\n", lo, hi,
              as.numeric(difftime(Sys.time(), t_mb, units = "mins")),
              hi, N_BOOT_METABRIC))
  utils::flush.console()
}

boot_metabric <- bind_rows(lapply(mb_files, readRDS)) %>% as_tibble()
saveRDS(boot_metabric, here::here("results/objects/bootstrap_metabric_prognostic.rds"))
cat("[TIMING] METABRIC bootstrap:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minutes\n")

# ------------------------------------------------------------------
# STEP 6. Equivalence testing + multiplicity
# ------------------------------------------------------------------
ep_map <- c(cindex_adjusted = "d_cindex_adj", cindex_cluster_only = "d_cindex_cl",
            lrt_over_age_er = "d_lrt")
sesoi_map <- c(cindex_adjusted = SESOI_CINDEX, cindex_cluster_only = SESOI_CINDEX,
               lrt_over_age_er = 3.84)

# Nothing below may abort: the bootstrap is finished and saved, and a
# formatting or plotting error must not destroy it.
equiv_metabric <- tryCatch(map_dfr(names(ep_map), function(ep) {
  col <- ep_map[[ep]]
  boot_metabric %>%
    group_by(k, view) %>%
    group_modify(~ {
      res <- tost_equivalence_bootstrap(.x[[col]], sesoi = sesoi_map[[ep]],
                                        label = paste(ep, unique(.y$view), sep = ":"))
      res$mdd_80pct_power <- minimum_detectable_difference(.x[[col]])
      res
    }) %>% ungroup() %>% mutate(endpoint = ep)
}) %>%
  group_by(k, endpoint) %>%
  group_modify(~ adjust_grid_multiplicity(.x)) %>%
  ungroup(),
  error = function(e) { cat("*** equivalence table FAILED:", conditionMessage(e),
    "\n    Raw bootstrap is safe in results/objects/ -- rebuildable in seconds.\n"); NULL })

if (!is.null(equiv_metabric)) {
  try(write_csv(equiv_metabric,
      here::here("results/tables/table_metabric_prognostic_added_value.csv")), silent = TRUE)
}

cat("\n\n=== HEADLINE: prognostic added value in METABRIC (adjusted C-index) ===\n")
if (!is.null(equiv_metabric)) try(print(equiv_metabric %>% filter(endpoint == "cindex_adjusted") %>%
        dplyr::select(k, view, point_estimate, ci95_lower, ci95_upper,
                      p_two_sided_BH, p_tost_BH, superior_after_BH,
                      equivalent_after_BH, mdd_80pct_power) %>%
        arrange(k, desc(point_estimate)), n = 40), silent = TRUE)

verdict_metabric <- tryCatch(equiv_metabric %>%
  group_by(endpoint, k) %>%
  summarise(
    n_views = n(),
    n_superior_after_BH = sum(superior_after_BH, na.rm = TRUE),
    n_equivalent_after_BH = sum(equivalent_after_BH, na.rm = TRUE),
    best_view = view[which.max(point_estimate)],
    best_point_estimate = max(point_estimate, na.rm = TRUE),
    .groups = "drop"
  ), error = function(e) { cat("*** verdict FAILED:", conditionMessage(e), "\n"); NULL })

if (!is.null(verdict_metabric)) {
  cat("\n=== METABRIC grid verdict ===\n")
  try(print(verdict_metabric, width = Inf), silent = TRUE)
  try(write_csv(verdict_metabric,
      here::here("results/tables/table_metabric_prognostic_verdict.csv")), silent = TRUE)
}

# Forest plot
p <- (if (is.null(equiv_metabric)) tibble() else equiv_metabric) %>%
  filter(endpoint == "cindex_adjusted") %>%
  ggplot(aes(x = reorder(view, point_estimate), y = point_estimate)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_hline(yintercept = c(-SESOI_CINDEX, SESOI_CINDEX), linetype = "dashed", colour = "firebrick") +
  geom_pointrange(aes(ymin = ci95_lower, ymax = ci95_upper)) +
  coord_flip() + facet_wrap(~ paste0("k = ", k)) +
  labs(x = NULL, y = "Delta C-index (fused - expression only), METABRIC",
       title = "Prognostic added value of TF-activity integration, external cohort",
       subtitle = "Adjusted for age + ER. Dashed red = pre-specified SESOI (0.02 C-index)") +
  theme_bw(base_size = 10)
try(ggsave(here::here("results/figures/figure_metabric_prognostic_added_value.png"),
           p, width = 11, height = 6, dpi = 150), silent = TRUE)

writeLines(c(
  "WHY THIS ENDPOINT, AND HOW TO REPORT IT.",
  "======================================================================",
  "Script 07 measures added value as agreement with PAM50. PAM50 labels",
  "are derived FROM gene expression, so an expression-only clustering is",
  "near-optimal for that target by construction and the baseline arm has",
  "a structural advantage the TF arm does not. A null result on such an",
  "endpoint is close to uninformative.",
  "",
  "This script measures added value as PROGNOSTIC DISCRIMINATION OF",
  "OVERALL SURVIVAL IN AN INDEPENDENT COHORT. Survival is not computed",
  "from the expression matrix, so neither arm is favoured; METABRIC took",
  "no part in deriving the clusters; and METABRIC supplies roughly an",
  "order of magnitude more events than TCGA with far longer follow-up.",
  "This is the primary added-value analysis for the manuscript. The",
  "PAM50 result should be reported as a secondary, explicitly circular",
  "comparison.",
  "",
  "SESOI JUSTIFICATION (C-index, 0.02).",
  "A 0.02 gain in C-index is at the lower edge of what is generally",
  "regarded as a clinically meaningful improvement in discrimination for",
  "a prognostic model. Choosing a small bound makes the equivalence test",
  "HARDER to pass, so passing it is a strong claim; and any gain below",
  "it could not change clinical decision-making regardless of its",
  "p-value.",
  "",
  "THE lrt_over_age_er ENDPOINT MATTERS INDEPENDENTLY.",
  "It measures information the clustering adds BEYOND age and ER status.",
  "A cluster solution that merely re-encodes ER scores near zero here.",
  "This is precisely the discrimination the TCGA-only survival analysis",
  "could not make, and it is the correct place to check whether the",
  "strong unadjusted METABRIC log-rank result reflects genuine",
  "subtype information or simply recapitulates ER status.",
  "",
  "LIMITATION -- BOOTSTRAP SCOPE.",
  "The bootstrap resamples METABRIC patients while holding the",
  "TCGA-derived centroids fixed. This is the correct frame for",
  "evaluating a fixed, pre-trained classifier on new data, and it is",
  "what a clinical reader assumes. It does not propagate uncertainty in",
  "the discovery-cohort clustering itself. A nested bootstrap that",
  "re-derived clusters within each replicate would, at a very large",
  "computational cost, and would widen these intervals.",
  "",
  "LIMITATION -- COVARIATE AVAILABILITY.",
  "This METABRIC export has no stage or grade column, so adjustment is",
  "age + ER only, whereas the TCGA discovery models adjust for age +",
  "stage. The adjustment sets are therefore NOT matched across cohorts.",
  "",
  "LIMITATION -- FEATURE SPACE.",
  "Projection uses the TCGA top-5000-variance gene set, which is a",
  "TCGA-optimised feature space, not an independently derived signature.",
  "Both arms inherit this equally, so it does not bias the CONTRAST,",
  "but it does limit the absolute performance of every arm."
), here::here("results/tables/NOTE_prognostic_added_value_design.txt"))

cat("\n\n================ PRODUCED OUTPUTS ================\n")
for (f in c("results/objects/bootstrap_metabric_prognostic.rds",
            "results/tables/table_metabric_prognostic_point_scores.csv",
            "results/tables/table_metabric_prognostic_added_value.csv",
            "results/tables/table_metabric_prognostic_verdict.csv",
            "results/figures/figure_metabric_prognostic_added_value.png")) {
  fp <- here::here(f)
  cat(sprintf("  %-8s %-56s %s\n", if (file.exists(fp)) "PRESENT" else "MISSING", f,
              if (file.exists(fp)) paste0(round(file.size(fp)/1024, 1), " KB") else ""))
}
cat(sprintf("  METABRIC bootstrap batches on disk: %d\n",
            length(list.files(here::here("results/objects"), "^metabric_boot_"))))
cat("==================================================\n\n")

log_session_info(script_name, key_packages = c("survival", "SNFtool", "org.Hs.eg.db"))
cat("\n✓ 14_prognostic_added_value.R complete.\n")

close_logger(log_con, script_name)
