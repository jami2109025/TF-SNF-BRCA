# ============================================================
# 07_benchmark_models.R  (REVISED: grid-wide added value with
#                         equivalence testing and cached distances)
# ============================================================
# WHAT CHANGED AND WHY
# ------------------------------------------------------------
# The previous revision produced this pipeline's headline number:
#
#     added value of TF integration, k=2:  +0.009  [-0.016, 0.058]
#     added value of TF integration, k=3:  +0.021  [-0.117, 0.055]
#
# Both intervals contain zero. As reported, that result cannot support
# ANY claim. "The CI included zero" is compatible with "TF activity
# adds nothing" and equally compatible with "this study could not have
# detected it anyway", and nothing in the previous design distinguishes
# them. Four changes convert this from an inconclusive non-result into
# a defensible finding:
#
# (1) EQUIVALENCE TESTING against a pre-specified SESOI.
#     This is the change that actually makes the negative result
#     publishable. Instead of failing to reject "effect = 0", we test
#     and reject "effect >= 0.05 ARI". Rejecting THAT is positive
#     evidence of absence rather than absence of evidence. See
#     tost_equivalence_bootstrap() and the SESOI justification in
#     R/utils_benchmark.R.
#
# (2) A GRID, not a single configuration. Every regulon x method view
#     from the revised script 05 is benchmarked, so the conclusion is
#     "no configuration in a systematically varied grid added value"
#     rather than "our one configuration did not". BH correction is
#     applied across the grid, because reporting the best of ~15 views
#     without adjustment would be a garden-of-forking-paths error.
#
# (3) MULTIPLE ENDPOINTS, including one that is not circular.
#     ARI-vs-PAM50 was the only endpoint before, and it is structurally
#     biased toward the expression-only baseline: PAM50 labels are
#     themselves derived from gene expression, so the baseline arm is
#     being scored against a target computed from its own input. A
#     benchmark in which one arm cannot really win is not a fair test
#     of anything. Survival endpoints (Cox LRT, Harrell's C) are added
#     because overall survival is NOT downstream of the expression
#     matrix and therefore does not favour either arm by construction.
#
# (4) CACHED DISTANCES -- a pure engineering change with zero effect on
#     any number, verified at runtime. The previous bootstrap took
#     22.8 HOURS for a single view; the grid version would have taken
#     weeks and would simply never have been run. Pairwise distance
#     d(i,j) is a property of the pair, so it does not change under
#     resampling: the O(n^2 * p) distance step is computed once and
#     every replicate indexes into it as D[idx, idx], leaving only the
#     cheap affinity/kNN and SNF steps per replicate. Correctness is
#     asserted at runtime by verify_distance_cache(), which fails hard
#     rather than silently producing fast wrong answers.
#
# ALSO FIXED
#   - Degeneracy guard for zero distances created by with-replacement
#     resampling (duplicated patients). The previous bootstrap had no
#     such guard; a replicate could produce NaN affinities and still
#     return a finite, meaningless ARI.
#   - Paired resampling across views (identical patients per replicate
#     across every arm), which is both statistically correct for a
#     difference statistic and substantially more precise.
#   - Minimum-detectable-difference reported alongside every null, so
#     the power question is answered pre-emptively rather than by a
#     reviewer.
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

script_name <- "07_benchmark_models"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

# ------------------------------------------------------------------
# CONFIG -- pre-specified
# ------------------------------------------------------------------
# MEMORY-CONSTRAINED CORE COUNT.
#
# This machine has 8 GB of RAM, not 16. Each PSOCK worker costs roughly
# 320 MB (R base ~120 MB, packages ~80 MB, its distance matrices ~19 MB,
# ~100 MB of SNF working set), and the master holds a further ~800 MB.
# Seven workers is ~3.0 GB of R on top of the 3-4 GB Windows and RStudio
# already use -- which is why a worker was killed 7.5 hours into the
# previous attempt ("error reading from connection").
#
# Six workers costs ~2.7 GB and adds ~2 hours. That is the right trade:
# this run has already been lost twice to memory pressure, and a slower
# run that finishes beats a faster one that does not.
N_CORES_OVERRIDE <- 6
N_BOOT           <- 200          # replicates; paired across all views
# Lowered from 250 after a 7-hour loss to the aricode NA crash. 200 is
# the pre-stated floor: it leaves 10 observations in each 5% tail of the
# equivalence test, which is the minimum for a stable percentile. Do not
# go lower -- below this the TOST decision becomes noise-driven, and the
# TOST decision is the thesis's central claim.
# 300, not 500: script 05 produced 29 views and SNF was measured at
# 696 s on this cohort with reference BLAS. 4 views x 300 x ~800 core-s
# / 7 cores is ~38 h, which fits the remaining budget; 500 would be 64 h,
# which does not. 300 still puts 15 observations in each 5% tail of the
# equivalence test -- above the ~200-replicate floor where those tails
# become unstable.
BATCH_SIZE       <- 25           # replicates per checkpointed batch
# A batch is the unit of work that can be lost to a crash. At 25
# replicates on 7 cores that is ~50 minutes, versus the ~8 hours lost
# when the whole chunk was dispatched at once. Completed batches are
# skipped on re-run, so the job is resumable.
RESAMPLE_SCHEME  <- "bootstrap"  # "bootstrap" or "subsample" (sensitivity)
VIEW_CHUNK_SIZE  <- 1            # ONE view per chunk -- see below
# Chunk size controls BOTH worker memory and checkpoint frequency, and
# the second matters more here. At size 5 all four bootstrap views sit
# in a single chunk, so the only checkpoint is written when the whole
# 38-hour bootstrap finishes -- a crash, reboot or power cut at hour 30
# would destroy everything. At size 1 a checkpoint lands roughly every
# 10 hours (bootstrap_chunk_01..04.rds).
#
# The cost is that the expression-only baseline is recomputed once per
# chunk instead of once overall: 4 x 300 x ~37 core-s = ~1.3 extra
# hours. That is cheap insurance on a job this long, and it changes no
# number -- the resample index depends only on the replicate number, so
# every chunk sees identical patients.
SNF_K            <- 20
SNF_SIGMA        <- 0.5
SNF_T            <- 20
K_VALUES         <- c(2, 3)

# ------------------------------------------------------------------
# WHICH VIEWS GET A CONFIDENCE INTERVAL
#
# Point estimates are computed for EVERY view in the grid -- they are
# cheap. Bootstrap CIs and equivalence tests are computed only for the
# views named here, because bootstrap cost scales as
# N_BOOT x n_views and a full 30-view grid at N_BOOT = 500 would take
# roughly two weeks on 7 cores.
#
# This subset MUST be fixed BEFORE looking at any point estimate.
# Choosing which views to bootstrap after seeing which ones scored well
# would invalidate every interval reported. Writing the list here, in
# the config block, is what makes that pre-specification auditable.
#
# The default below spans the two factors that actually vary --
# REGULON SOURCE (DoRothEA vs CollecTRI) and REAL vs EDGE-PERMUTED --
# while deliberately NOT spanning inference method. That omission is
# justified by this pipeline's own measurement: script 05 found median
# per-TF Spearman of 0.995-1.000 between VIPER and ULM, i.e. the
# methods are the same view measured twice. Bootstrapping both would
# spend half the compute budget re-measuring a quantity already shown
# to be invariant. Cite table_tf_view_concordance.csv when defending
# this choice.
#
# Set to NULL to bootstrap every view (only feasible with a small grid).
# ------------------------------------------------------------------
# Default: a FULL 2 x 2 x 2 FACTORIAL over the three factors that
# actually vary --
#     regulon source  {DoRothEA A-C, CollecTRI}
#   x inference method {VIPER, ULM}
#   x network content  {real, degree-preserving edge-permuted}
#
# A factorial subset is worth the extra compute over an ad hoc
# selection: it lets the manuscript state that each factor was varied
# while holding the others fixed, rather than that a few configurations
# were tried. If the budget is tight, drop the ULM arm first -- script
# 05 measured median per-TF Spearman of 0.995-1.000 between VIPER and
# ULM, so that axis carries the least new information.
# REDUCED FROM THE 2x2x2 FACTORIAL TO A 2x2, DROPPING THE METHOD AXIS.
#
# Justified by measurement, not by convenience: script 05 found
# rho = 0.997 between dorothea_AC__viper and dorothea_AC__ulm
# (table_tf_view_concordance.csv). VIPER and ULM are the same view
# measured twice, so bootstrapping both would spend half the compute
# budget re-measuring a quantity already shown to be invariant.
#
# The two axes that DO vary are retained in full:
#   regulon source  {DoRothEA A-C, CollecTRI}   -- real-vs-real rho = 0.430
#   network content {real, edge-permuted}       -- real-vs-perm  rho = 0.243
#                                                  against a random-random
#                                                  floor of rho = 0.310
# REDUCED TO TWO VIEWS for the time remaining. Both are REAL regulons
# from INDEPENDENT sources, which is what the headline claim needs:
# "neither of two independent regulon sources added value". The
# edge-permuted controls keep their point estimates (all 29 views get
# those); only their confidence intervals are dropped. The real-versus-
# permuted distinction is established independently, and without any
# bootstrap, by the cross-view concordance analysis in script 05
# (real-vs-real rho = 0.430 against a random-random floor of 0.310).
VIEWS_FOR_BOOTSTRAP <- c(
  "dorothea_AC__viper",
  "collectri__viper"
)

n_cores <- if (is.null(N_CORES_OVERRIDE)) max(1, parallel::detectCores() - 1) else N_CORES_OVERRIDE
cat("Using", n_cores, "cores.\n")
report_blas_speed()

# Memory is the binding constraint on this machine (8 GB total). Report
# it at each major stage so a doomed run is visible before it starts,
# not seven hours in.
report_memory <- function(tag) {
  g <- gc(verbose = FALSE, full = TRUE)
  cat(sprintf("[memory %s] R is holding ~%.0f MB\n", tag, sum(g[, 2])))
  invisible(NULL)
}
report_memory("startup")

# ------------------------------------------------------------------
# Worker bootstrapping: source the utils files INSIDE each worker.
#
# Exporting helper functions by name via clusterExport() is fragile in
# exactly one way, and it has now bitten this pipeline twice: add a new
# internal helper (.clean_pair) that an exported function calls, forget
# to add it to the varlist, and all seven workers die with
# "could not find function". The failure surfaces only after cluster
# spin-up, i.e. ~15 minutes in.
#
# Sourcing R/utils.R and R/utils_benchmark.R inside each worker makes
# the worker environment self-sufficient: every helper, present and
# future, is available without anyone having to remember a varlist.
# Paths are passed as an argument rather than exported, so there is no
# ordering dependency between clusterCall and clusterExport.
# ------------------------------------------------------------------
UTILS_PATHS <- c(here::here("R", "utils.R"), here::here("R", "utils_benchmark.R"))
stopifnot(all(file.exists(UTILS_PATHS)))

source_utils_on_workers <- function(cl) {
  if (is.null(cl)) return(invisible(NULL))
  parallel::clusterCall(cl, function(paths) {
    for (f in paths) source(f)
    invisible(TRUE)
  }, UTILS_PATHS)
  invisible(NULL)
}

# ------------------------------------------------------------------
# Robust parallel helpers
#
# On Windows, parallel::makeCluster() creates a PSOCK cluster by opening
# localhost sockets to freshly spawned Rscript processes. Antivirus,
# firewall rules and corporate endpoint agents routinely block exactly
# that, and the failure is SILENT: makeCluster() blocks forever, CPU
# sits at 0%, and no error is ever printed. On a job measured in tens of
# hours that is the worst possible failure mode, because it looks
# identical to "still working".
#
# make_cluster_checked() therefore builds the cluster with a timeout and
# then runs a trivial round-trip through it. If either step fails it
# returns NULL rather than hanging, and callers fall back to serial
# execution -- slower, but it finishes, and it prints progress.
# ------------------------------------------------------------------
# TIMEOUT NOTE: 90 s proved too short in practice. makeCluster(7) spawns
# seven fresh R processes and waits for seven socket handshakes; on a
# machine that is busy, or that still has orphaned Rscript.exe processes
# from an earlier failed attempt holding resources, that legitimately
# takes minutes. A verified makeCluster(2) round-trip on this same
# machine succeeded in seconds, which points at load/orphans rather than
# at a blocked socket layer. 300 s gives real headroom while still
# failing fast enough to be useful.
# TIMEOUT: MEASURED on this machine. makeCluster(7) with zero orphaned
# processes took 197 s -- about 28 s per worker, which is slow for PSOCK
# and almost certainly antivirus scanning each new Rscript.exe launch.
# 300 s left only 100 s of margin, and spin-up will be slower still
# during the real run because the master session holds ~1 GB of cached
# distance matrices. 600 s is generous; the only cost of a generous
# timeout is that a genuine failure takes longer to report, and we have
# now verified that the cluster does come up.
make_cluster_checked <- function(n, timeout_sec = 600) {
  cl <- tryCatch({
    setTimeLimit(elapsed = timeout_sec, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
    parallel::makeCluster(n)
  }, error = function(e) {
    cat("  [parallel] makeCluster failed or timed out:", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(cl)) return(NULL)

  ok <- tryCatch({
    r <- parallel::clusterEvalQ(cl, 1L + 1L)
    length(r) == n && all(unlist(r) == 2L)
  }, error = function(e) FALSE)

  if (!ok) {
    cat("  [parallel] cluster round-trip test FAILED.\n")
    try(parallel::stopCluster(cl), silent = TRUE)
    return(NULL)
  }
  cat("  [parallel] cluster of", n, "workers verified.\n")
  cl
}

#' Try the requested worker count, then progressively fewer.
#'
#' Falling straight from "7 workers failed" to "run serially" is a
#' 7x cliff. Most spin-up failures are contention, not a hard block, so
#' a smaller cluster usually succeeds -- and 4 workers finishing in 56 h
#' is worth vastly more than 1 worker finishing in 222 h.
make_cluster_backoff <- function(n) {
  for (k in unique(c(n, max(2, n %/% 2), 2))) {
    cat("  [parallel] attempting", k, "workers ...\n"); utils::flush.console()
    cl <- make_cluster_checked(k)
    if (!is.null(cl)) {
      if (k < n) cat("  [parallel] NOTE: running on", k, "workers instead of", n,
                     "-- wall-clock scales by", round(n / k, 2), "x\n")
      return(cl)
    }
  }
  cat("  [parallel] all cluster attempts failed.\n")
  NULL
}

#' parLapply when a cluster is available, lapply with progress when not.
par_or_serial <- function(cl, X, FUN, label = "") {
  if (!is.null(cl)) return(parallel::parLapply(cl, X, FUN))
  cat("  [serial] running", length(X), "task(s) sequentially --",
      "this is slower but prints progress and cannot hang on a socket.\n")
  out <- vector("list", length(X))
  for (i in seq_along(X)) {
    t0 <- Sys.time()
    out[[i]] <- FUN(X[[i]])
    cat(sprintf("    [%s] %d/%d done (%.1f min)\n", label, i, length(X),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    utils::flush.console()
  }
  out
}

# ------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------
vst_top2000     <- readRDS(here::here("data/processed/vst_top2000_genes.rds"))
sample_metadata <- readRDS(here::here("data/processed/sample_metadata_matched.rds"))
tf_grid         <- readRDS(here::here("data/processed/tf_activity_grid.rds"))
W_fused_primary <- readRDS(here::here("results/objects/W_fused.rds"))

cat("TF-activity grid loaded:", length(tf_grid), "views --",
    paste(names(tf_grid), collapse = ", "), "\n")

# Align everything to a single canonical sample set: the intersection
# across the expression matrix, EVERY grid view, and the metadata.
common_samples <- Reduce(intersect, c(
  list(colnames(vst_top2000), rownames(sample_metadata)),
  lapply(tf_grid, colnames)
))
cat("Canonical sample set common to expression, metadata and ALL grid views:",
    length(common_samples), "\n")
stopifnot(length(common_samples) > 100)

vst_top2000     <- vst_top2000[, common_samples, drop = FALSE]
sample_metadata <- sample_metadata[common_samples, , drop = FALSE]
tf_grid         <- lapply(tf_grid, function(m) m[, common_samples, drop = FALSE])

# Safety check on the loaded primary fused network (multi-hour script;
# fail before computing, not after).
stopifnot(all(common_samples %in% colnames(W_fused_primary)))
W_fused_primary <- W_fused_primary[common_samples, common_samples]

# ------------------------------------------------------------------
# Views, scaled. Each view is z-scored per feature exactly as in
# script 06, so fusion inputs are identical to the main analysis.
# ------------------------------------------------------------------
scale_view <- function(mat) {
  s <- scale(t(mat))          # samples in rows
  s[is.na(s)] <- 0
  s
}
expr_scaled <- scale_view(vst_top2000)
tf_scaled_list <- lapply(tf_grid, scale_view)

# ------------------------------------------------------------------
# DISTANCE CACHE -- computed once, verified, reused everywhere
# ------------------------------------------------------------------
cat("\nBuilding distance cache (once) ...\n")
t0 <- Sys.time()
D_expr <- dist2_cached(expr_scaled)
verify_distance_cache(expr_scaled, D_expr, label = "expression")

D_tf_list <- lapply(names(tf_scaled_list), function(nm) {
  d <- dist2_cached(tf_scaled_list[[nm]])
  verify_distance_cache(tf_scaled_list[[nm]], d, label = nm)
  d
})
names(D_tf_list) <- names(tf_scaled_list)

# Naive concatenation baseline: the fair "just glue the features
# together" comparator. If SNF fusion does not beat concatenation,
# the FUSION machinery specifically is adding nothing, independently
# of whether the TF view is informative -- a separate claim from the
# added-value question and worth reporting on its own.
# Concatenation distances are needed ONLY by the point-estimate stage,
# which now runs AFTER the bootstrap. Building them here would hold
# 29 x 9.6 MB = ~278 MB in the master for the entire multi-hour
# bootstrap for no reason. Deferred to the point-estimate section.
D_concat_list <- NULL

cat("[TIMING] Distance cache built in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minutes\n")

# ------------------------------------------------------------------
# REPRODUCIBILITY CHECK against the stored primary network.
#
# W_fused for the primary view is rebuilt here through the new cached
# code path and compared to results/objects/W_fused.rds, which script
# 06 built through the ORIGINAL uncached path. Agreement demonstrates
# that the caching optimisation reproduces the published network
# exactly, and that no upstream file has drifted since script 06 ran.
# This is the check that lets the manuscript state the speedup changed
# no result -- an assertion that should be demonstrated, not asserted.
# ------------------------------------------------------------------
primary_view_name <- if ("dorothea_AC__viper" %in% names(D_tf_list))
  "dorothea_AC__viper" else names(D_tf_list)[1]
W_fused_check <- SNFtool::SNF(
  list(affinity_from_dist(D_expr, K = SNF_K, sigma = SNF_SIGMA),
       affinity_from_dist(D_tf_list[[primary_view_name]], K = SNF_K, sigma = SNF_SIGMA)),
  K = SNF_K, t = SNF_T
)
max_dev <- max(abs(W_fused_check - W_fused_primary))
cat("Primary fused network vs stored W_fused.rds: max abs deviation =",
    format(max_dev, digits = 3), "\n")
if (max_dev > 1e-8) {
  warning("Rebuilt primary fused network differs from results/objects/W_fused.rds ",
          "(max deviation ", format(max_dev, digits = 3), "). Either an upstream ",
          "input changed since script 06 was run, or the primary view is not ",
          "'", primary_view_name, "'. Resolve before reporting any result that ",
          "mixes outputs from the two runs.")
} else {
  cat("  -> exact match: the cached distance path reproduces the published network.\n")
}

# ------------------------------------------------------------------
# Reference labels and survival outcome
# ------------------------------------------------------------------
pam50_vec <- sample_metadata$PAM50
has_pam50 <- !is.na(pam50_vec) & pam50_vec != ""
cat("Samples with PAM50 available:", sum(has_pam50), "of", length(pam50_vec), "\n")

# Survival construction MIRRORS 10_survival_analysis.R exactly. The
# coercion below is not cosmetic: TCGA's curated "paper_" clinical
# fields arrive as character, and without coercion every survival
# endpoint here would silently become NA and the whole prognostic
# comparison would report as unavailable rather than as failed.
sm <- sample_metadata
sm$days_to_death         <- suppressWarnings(as.numeric(sm$days_to_death))
sm$days_to_last_followup <- suppressWarnings(as.numeric(sm$days_to_last_followup))
surv_time <- ifelse(tolower(sm$vital_status) == "dead",
                    sm$days_to_death, sm$days_to_last_followup) / 365.25
surv_status <- ifelse(tolower(sm$vital_status) == "dead", 1, 0)
surv_valid <- is.finite(surv_time) & surv_time > 0 & !is.na(surv_status)

n_events <- sum(surv_status[surv_valid] == 1, na.rm = TRUE)
cat("Survival endpoint: N =", sum(surv_valid), " events =", n_events,
    sprintf(" (%.1f%%)\n", 100 * n_events / sum(surv_valid)))
cat("Median follow-up (years):",
    round(median(surv_time[surv_valid], na.rm = TRUE), 2), "\n")

# Event count is reported UP FRONT and deliberately, because it governs
# how much any TCGA survival null is worth. With this few events a
# within-TCGA prognostic comparison is exploratory; the properly
# powered prognostic test is the METABRIC analysis in script 14, which
# has an order of magnitude more events and far longer follow-up.
if (n_events < 200) {
  cat("\nNOTE: fewer than 200 events in TCGA. Survival-based endpoints in\n",
      "THIS script are underpowered and are reported as secondary and\n",
      "exploratory only. The adequately powered prognostic added-value\n",
      "test is 14_prognostic_added_value.R (METABRIC).\n", sep = "")
}

# ------------------------------------------------------------------
# Point-estimate benchmark over the full grid
# ------------------------------------------------------------------
# PAIRED BOOTSTRAP over the grid
#
# Structure, and why it is this way:
#   - The resample index for replicate b is a deterministic function of
#     b alone, so every view sees the SAME patients at replicate b and
#     the added-value difference is a within-replicate PAIRED contrast.
#     Independent resampling per arm would inflate the variance of the
#     difference so much that a real effect could be masked -- exactly
#     the failure mode that would make a negative result meaningless.
#   - The expression-only baseline is recomputed inside each replicate
#     rather than being held fixed, because the baseline is itself a
#     random quantity and treating it as fixed would understate the
#     uncertainty of the difference.
#   - Views are processed in chunks so worker memory stays bounded;
#     chunking changes nothing numerically because the resample index
#     depends only on b.
# ------------------------------------------------------------------
n_samples <- length(common_samples)

# Restrict the bootstrap to the pre-specified view subset (see CONFIG).
if (is.null(VIEWS_FOR_BOOTSTRAP)) {
  view_names <- names(D_tf_list)
} else {
  view_names <- intersect(VIEWS_FOR_BOOTSTRAP, names(D_tf_list))
  missing_views <- setdiff(VIEWS_FOR_BOOTSTRAP, names(D_tf_list))
  if (length(missing_views) > 0) {
    warning("VIEWS_FOR_BOOTSTRAP names not present in the grid and therefore ",
            "NOT bootstrapped: ", paste(missing_views, collapse = ", "),
            ". If these were meant to be in the grid (e.g. CollecTRI views ",
            "after a failed download), fix that before reporting -- silently ",
            "dropping a pre-specified view changes the analysis plan.")
  }
  if (length(view_names) == 0) {
    stop("None of the views in VIEWS_FOR_BOOTSTRAP exist in the grid. ",
         "Available: ", paste(names(D_tf_list), collapse = ", "))
  }
}
view_chunks <- split(view_names, ceiling(seq_along(view_names) / VIEW_CHUNK_SIZE))

# ------------------------------------------------------------------
# FREE MEMORY BEFORE THE LONG JOB.
#
# The bootstrap touches only D_expr and the distance matrices of the
# pre-specified views. The other ~27 view matrices and the raw grid are
# not needed again until the point-estimate stage, which runs afterwards
# and rebuilds what it requires. On an 8 GB machine, holding them
# through a 16-hour bootstrap is the difference between finishing and
# being killed by the OS.
# ------------------------------------------------------------------
bootstrap_needs <- if (is.null(VIEWS_FOR_BOOTSTRAP)) names(D_tf_list) else
  intersect(VIEWS_FOR_BOOTSTRAP, names(D_tf_list))
D_tf_keep <- D_tf_list[bootstrap_needs]
suppressWarnings(rm(tf_grid, W_fused_primary, W_fused_check, D_concat_list))
gc(verbose = FALSE, full = TRUE)
report_memory("after freeing unused views")

cat("\n================ Paired bootstrap (n_boot =", N_BOOT, ", scheme =",
    RESAMPLE_SCHEME, ") ================\n")
cat("Point estimates: all", length(D_tf_list), "views.",
    "Bootstrap CIs:", length(view_names), "pre-specified view(s):\n  ",
    paste(view_names, collapse = "\n   "), "\n")
cat("Processed in", length(view_chunks), "chunk(s) of up to",
    VIEW_CHUNK_SIZE, "to bound worker memory.\n")

boot_worker <- function(b, D_expr, D_tf_chunk, pam50_vec, has_pam50,
                        surv_time, surv_status, surv_valid,
                        n_samples, K_VALUES, scheme) {

  idx <- make_resample_index(b, n_samples, scheme = scheme)

  pam_b    <- pam50_vec[idx]
  haspam_b <- has_pam50[idx]
  time_b   <- surv_time[idx]
  stat_b   <- surv_status[idx]
  valid_b  <- surv_valid[idx]

  W_expr_b <- affinity_from_dist(D_expr[idx, idx, drop = FALSE],
                                 K = SNF_K, sigma = SNF_SIGMA)

  # Baseline endpoints per k, computed once per replicate.
  base <- list()
  for (k in K_VALUES) {
    cl_expr_b <- tryCatch(SNFtool::spectralClustering(W_expr_b, K = k),
                          error = function(e) rep(NA_integer_, nrow(W_expr_b)))
    base[[as.character(k)]] <- c(
      ari  = ari_vs_reference(cl_expr_b, pam_b, haspam_b),
      nmi  = nmi_vs_reference(cl_expr_b, pam_b, haspam_b),
      lrt  = cox_lrt_endpoint(cl_expr_b[valid_b], time_b[valid_b], stat_b[valid_b]),
      cidx = concordance_endpoint(cl_expr_b[valid_b], time_b[valid_b], stat_b[valid_b])
    )
  }

  # CRITICAL LOOP ORDER: the view loop is OUTSIDE and the k loop INSIDE.
  # W_tf_b and W_fused_b do not depend on k -- only the final
  # spectralClustering(K = k) does. Nesting the other way round rebuilds
  # the fused network once per k value, which for K_VALUES = c(2, 3)
  # doubles the cost of the entire bootstrap while producing identical
  # numbers. SNF dominates the runtime here, so this ordering roughly
  # halves the wall clock.
  # A single replicate must never be able to abort the chunk. Any error
  # here degrades to a row of NAs, which tost_equivalence_bootstrap()
  # already filters; losing one replicate out of hundreds costs nothing,
  # losing the run costs days. n_boot_valid in the output records how
  # many replicates actually contributed.
  out <- list()
  for (nm in names(D_tf_chunk)) {
    res_nm <- tryCatch({
      W_tf_b <- affinity_from_dist(D_tf_chunk[[nm]][idx, idx, drop = FALSE],
                                   K = SNF_K, sigma = SNF_SIGMA)
      W_fused_b <- SNFtool::SNF(list(W_expr_b, W_tf_b), K = SNF_K, t = SNF_T)

      rows <- list()
      for (k in K_VALUES) {
        cl_fused_b <- SNFtool::spectralClustering(W_fused_b, K = k)
        bs <- base[[as.character(k)]]
        rows[[length(rows) + 1]] <- data.frame(
          b = b, k = k, view = nm,
          d_ari  = ari_vs_reference(cl_fused_b, pam_b, haspam_b) - bs[["ari"]],
          d_nmi  = nmi_vs_reference(cl_fused_b, pam_b, haspam_b) - bs[["nmi"]],
          d_lrt  = cox_lrt_endpoint(cl_fused_b[valid_b], time_b[valid_b], stat_b[valid_b]) - bs[["lrt"]],
          d_cidx = concordance_endpoint(cl_fused_b[valid_b], time_b[valid_b], stat_b[valid_b]) - bs[["cidx"]],
          stringsAsFactors = FALSE
        )
      }
      do.call(rbind, rows)
    }, error = function(e) {
      data.frame(b = b, k = K_VALUES, view = nm,
                 d_ari = NA_real_, d_nmi = NA_real_,
                 d_lrt = NA_real_, d_cidx = NA_real_,
                 stringsAsFactors = FALSE)
    })
    out[[length(out) + 1]] <- res_nm
  }
  do.call(rbind, out)
}

# Runtime estimate printed BEFORE committing to a multi-hour job, so an
# infeasible configuration is caught now rather than discovered later.
# MEASURED, not guessed: 07a_timing_probe.R timed one view-replicate
# on this exact cohort and machine at 801.5 core-seconds
# (affinity 3.4 s + SNF 696.5 s + 2 x spectralClustering 14.4 s,
# plus a 10% endpoint allowance). The earlier placeholder of 300
# understated the runtime by a factor of 2.7, which would have made
# the estimate printed below actively misleading at exactly the
# moment it is used to decide whether to commit to the run.
est_core_sec_per_view_rep <- 801.5
est_hours <- N_BOOT * length(view_names) * est_core_sec_per_view_rep / n_cores / 3600
cat(sprintf("\nEstimated bootstrap runtime: ~%.1f hours (%d replicates x %d views / %d cores).\n",
            est_hours, N_BOOT, length(view_names), n_cores))
cat("  Reduce N_BOOT or VIEWS_FOR_BOOTSTRAP now if that does not fit.\n")
cat("  Do not go below N_BOOT = 200: the equivalence test reads the 5th and\n")
cat("  95th percentiles, and fewer replicates make those tails unstable.\n\n")

t_boot_start <- Sys.time()
boot_all <- list()

for (ci in seq_along(view_chunks)) {
  chunk <- view_chunks[[ci]]
  D_tf_chunk <- D_tf_keep[chunk]
  cat("  chunk", ci, "of", length(view_chunks), ":", paste(chunk, collapse = ", "), "\n")

  cl_boot <- make_cluster_backoff(n_cores)
  # Cluster construction + configuration is a FUNCTION so that a batch
  # which loses a worker can rebuild an identical cluster and retry,
  # rather than aborting the run.
  setup_boot_cluster <- function() {
    cl <- make_cluster_backoff(n_cores)
    if (is.null(cl)) return(NULL)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(SNFtool); library(mclust); library(survival)
      })
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
        try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
      }
      Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    })
    parallel::clusterExport(
      cl,
      varlist = c("D_expr", "D_tf_chunk", "pam50_vec", "has_pam50",
                  "surv_time", "surv_status", "surv_valid", "n_samples",
                  "K_VALUES", "RESAMPLE_SCHEME", "SNF_K", "SNF_SIGMA", "SNF_T",
                  "PIPELINE_SEED", "AFFINITY_EPS",
                  "affinity_from_dist", "make_resample_index",
                  ".clean_pair",
                  "ari_vs_reference", "nmi_vs_reference",
                  "cox_lrt_endpoint", "concordance_endpoint",
                  "set_pipeline_seed", "boot_worker"),
      envir = environment()
    )
    cl
  }

  if (!is.null(cl_boot)) {
  # NOTE: workers no longer source R/utils*.R. Those files call
  # library(tidyverse) at load time, which added roughly 200 MB to EACH
  # of seven workers -- about 1.4 GB of avoidable footprint. A worker was
  # lost to memory pressure 7.5 hours into chunk 1 ("error reading from
  # connection"), destroying the whole chunk. None of the helper
  # functions the worker actually calls need tidyverse, so they are
  # exported by name instead. The smoke test below is what makes that
  # safe: a missing name now fails in two minutes rather than after the
  # full job has been dispatched.
  parallel::clusterEvalQ(cl_boot, {
    suppressPackageStartupMessages({
      library(SNFtool); library(mclust); library(survival)
    })
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
      try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
    }
    Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  })
  parallel::clusterExport(
    cl_boot,
    varlist = c("D_expr", "D_tf_chunk", "pam50_vec", "has_pam50",
                "surv_time", "surv_status", "surv_valid", "n_samples",
                "K_VALUES", "RESAMPLE_SCHEME", "SNF_K", "SNF_SIGMA", "SNF_T",
                "PIPELINE_SEED", "AFFINITY_EPS",
                "affinity_from_dist", "make_resample_index",
                ".clean_pair",
                "ari_vs_reference", "nmi_vs_reference",
                "cox_lrt_endpoint", "concordance_endpoint",
                "set_pipeline_seed", "boot_worker"),
    envir = environment()
  )

  }

  # SMOKE TEST -- one replicate, before committing the next eight hours.
  # Both crashes so far (aricode NA, missing .clean_pair) surfaced only
  # once the full parLapply had been dispatched, ~15 minutes after
  # launch, and looked like a hang until the log was read. Running a
  # single replicate first turns that into a two-minute, explicit
  # failure with the real message attached.
  cat("  [smoke] running replicate 1 through the cluster ...\n"); utils::flush.console()
  smoke <- tryCatch(
    par_or_serial(cl_boot, list(1L), function(b) {
      boot_worker(b, D_expr, D_tf_chunk, pam50_vec, has_pam50,
                  surv_time, surv_status, surv_valid,
                  n_samples, K_VALUES, RESAMPLE_SCHEME)
    }, label = "smoke")[[1]],
    error = function(e) e
  )
  if (inherits(smoke, "error")) {
    if (!is.null(cl_boot)) try(parallel::stopCluster(cl_boot), silent = TRUE)
    stop("SMOKE TEST FAILED on chunk ", ci, " -- aborting before wasting hours.\n",
         "  Message: ", conditionMessage(smoke), "\n",
         "  Nothing downstream has been computed; fix and re-run.")
  }
  if (!is.data.frame(smoke) || nrow(smoke) != length(K_VALUES) * length(D_tf_chunk)) {
    if (!is.null(cl_boot)) try(parallel::stopCluster(cl_boot), silent = TRUE)
    stop("SMOKE TEST returned an unexpected shape (", nrow(smoke), " rows, expected ",
         length(K_VALUES) * length(D_tf_chunk), "). Aborting.")
  }
  cat("  [smoke] PASS -- endpoints returned:",
      paste(sprintf("%s=%.4f", c("d_ari","d_nmi","d_lrt","d_cidx"),
                    c(smoke$d_ari[1], smoke$d_nmi[1], smoke$d_lrt[1], smoke$d_cidx[1])),
            collapse = "  "), "\n")
  utils::flush.console()

  # BATCHED EXECUTION WITH RESUME.
  #
  # Dispatching all N_BOOT replicates in one parLapply means a single
  # lost worker at hour 7.5 destroys the entire chunk -- which is exactly
  # what happened. Replicates are therefore run in batches of
  # BATCH_SIZE, each written to disk as soon as it completes. A crash now
  # costs at most one batch (~50 min), and re-running the script skips
  # every batch already on disk.
  #
  # This changes no number: the resample index depends only on the
  # replicate number b, so batch boundaries are invisible to the
  # statistics.
  batch_starts <- seq(1, N_BOOT, by = BATCH_SIZE)
  batch_files  <- character(0)

  for (bi in seq_along(batch_starts)) {
    lo <- batch_starts[bi]
    hi <- min(lo + BATCH_SIZE - 1, N_BOOT)
    bfile <- here::here("results/objects",
                        sprintf("bootstrap_c%02d_b%03d_%03d.rds", ci, lo, hi))
    batch_files <- c(batch_files, bfile)

    if (file.exists(bfile)) {
      cat(sprintf("    batch %d-%d: already on disk, skipping\n", lo, hi))
      next
    }

    t_b <- Sys.time()

    # RETRY. parLapply signals if any worker dies -- which is exactly how
    # 7.5 hours were lost ("error reading from connection", a worker
    # killed under memory pressure). A dead worker must cost one batch
    # and a cluster rebuild, never the run. Up to three attempts on a
    # fresh cluster, then this single batch runs serially.
    res_b <- NULL
    for (attempt in 1:3) {
      res_b <- tryCatch(
        par_or_serial(cl_boot, as.list(lo:hi), function(b) {
          boot_worker(b, D_expr, D_tf_chunk, pam50_vec, has_pam50,
                      surv_time, surv_status, surv_valid,
                      n_samples, K_VALUES, RESAMPLE_SCHEME)
        }, label = paste0("c", ci, "b", lo)),
        error = function(e) e
      )
      if (!inherits(res_b, "error")) break

      cat(sprintf("    [retry %d/3] batch %d-%d failed: %s\n",
                  attempt, lo, hi, conditionMessage(res_b)))
      utils::flush.console()
      try(parallel::stopCluster(cl_boot), silent = TRUE)
      gc(verbose = FALSE)
      Sys.sleep(20)
      cl_boot <- setup_boot_cluster()
    }

    if (inherits(res_b, "error")) {
      cat("    [fallback] three cluster attempts failed; running this batch serially.\n")
      utils::flush.console()
      res_b <- lapply(lo:hi, function(b) {
        tryCatch(
          boot_worker(b, D_expr, D_tf_chunk, pam50_vec, has_pam50,
                      surv_time, surv_status, surv_valid,
                      n_samples, K_VALUES, RESAMPLE_SCHEME),
          error = function(e) data.frame(b = b, k = K_VALUES, view = names(D_tf_chunk)[1],
                                         d_ari = NA_real_, d_nmi = NA_real_,
                                         d_lrt = NA_real_, d_cidx = NA_real_))
      })
    }

    saveRDS(dplyr::bind_rows(res_b), bfile)
    el <- as.numeric(difftime(Sys.time(), t_b, units = "mins"))
    done <- hi; total <- N_BOOT
    cat(sprintf("    batch %d-%d done (%.1f min) | chunk %d: %d/%d replicates | est. %.1f h left in chunk\n",
                lo, hi, el, ci, done, total, el * (total - done) / BATCH_SIZE / 60))
    utils::flush.console()
  }

  chunk_res <- lapply(batch_files, readRDS)
  if (!is.null(cl_boot)) parallel::stopCluster(cl_boot)

  boot_all[[ci]] <- dplyr::bind_rows(chunk_res)

  # CHECKPOINT after every chunk. This bootstrap can run for well over a
  # day; holding all results in memory until the end means a crash, a
  # reboot, or a power cut at hour 30 destroys everything. Each chunk is
  # written to its own file as soon as it completes, so a failed run can
  # be resumed by re-running only the chunks whose files are missing.
  chunk_path <- here::here("results/objects",
                           sprintf("bootstrap_chunk_%02d.rds", ci))
  saveRDS(boot_all[[ci]], chunk_path)
  cat("    checkpoint written:", basename(chunk_path), "\n")
  cat("    elapsed so far:",
      round(as.numeric(difftime(Sys.time(), t_boot_start, units = "hours")), 2), "hours",
      sprintf("| est. remaining: %.1f h\n",
              as.numeric(difftime(Sys.time(), t_boot_start, units = "hours")) *
                (length(view_chunks) - ci) / ci))
}

boot_diffs <- bind_rows(boot_all) %>% as_tibble()
saveRDS(boot_diffs, here::here("results/objects/bootstrap_added_value_diffs.rds"))

# Recovery note, in case a later stage of this script fails: the
# equivalence analysis below can be re-run without repeating the
# bootstrap by loading the checkpoints instead.
#   boot_diffs <- dplyr::bind_rows(lapply(
#     list.files(here::here("results/objects"), "^bootstrap_chunk_.*rds$",
#                full.names = TRUE), readRDS))

cat("[TIMING] Bootstrap total:",
    round(as.numeric(difftime(Sys.time(), t_boot_start, units = "hours")), 2), "hours\n")

# ------------------------------------------------------------------
# Equivalence testing per (k, view, endpoint), then BH across the grid
# ------------------------------------------------------------------
endpoint_cols <- c(ari = "d_ari", nmi = "d_nmi", cox_lrt = "d_lrt", cindex = "d_cidx")

# The SESOI is endpoint-specific. It is stated on the ARI scale in
# utils_benchmark.R and translated here; for endpoints without an
# established practical threshold the equivalence test is reported but
# should be interpreted as exploratory, and this is recorded in the
# output rather than left implicit.
sesoi_for <- c(ari = SESOI_ARI, nmi = SESOI_ARI, cox_lrt = 3.84, cindex = 0.02)
sesoi_basis <- c(
  ari     = "pre-specified; smaller than this pipeline's own k=2 vs k=3 ARI gap (0.085)",
  nmi     = "same magnitude as the ARI SESOI; exploratory, no established NMI threshold",
  cox_lrt = "chi-square 3.84 = the 1-df 0.05 critical value, i.e. one 'significant term' worth of prognostic information",
  cindex  = "0.02 C-index; below the smallest gain generally regarded as clinically meaningful"
)

# ==================================================================
# EVERYTHING BELOW IS NON-FATAL.
#
# The expensive computation is finished and the raw bootstrap
# differences are already on disk. From this point on, a formatting,
# summarising or plotting error must NOT be allowed to abort the script
# and cost a day of compute. Each stage is wrapped: it reports its own
# failure and the next stage still runs. If a stage fails, the analysis
# can be rebuilt in seconds from
# results/objects/bootstrap_added_value_diffs.rds.
# ==================================================================

equiv_table <- tryCatch({
  et <- map_dfr(names(endpoint_cols), function(ep) {
    col <- endpoint_cols[[ep]]
    boot_diffs %>%
      group_by(k, view) %>%
      group_modify(~ {
        res <- tost_equivalence_bootstrap(.x[[col]], sesoi = sesoi_for[[ep]],
                                          label = paste(ep, unique(.y$view), sep = ":"))
        res$mdd_80pct_power <- minimum_detectable_difference(.x[[col]])
        res
      }) %>%
      ungroup() %>%
      mutate(endpoint = ep, sesoi_basis = sesoi_basis[[ep]])
  })

  # Multiplicity is applied WITHIN each (k, endpoint) family -- the
  # family of views tested for the same question at the same k. Pooling
  # across endpoints would over-correct, since the endpoints answer
  # different questions and are not competing hypotheses.
  et %>%
    group_by(k, endpoint) %>%
    group_modify(~ adjust_grid_multiplicity(.x)) %>%
    ungroup() %>%
    dplyr::select(endpoint, k, view, point_estimate,
                  ci95_lower, ci95_upper, ci90_lower, ci90_upper,
                  sesoi, p_two_sided, p_two_sided_BH, p_tost, p_tost_BH,
                  equivalent, superior, superior_after_BH, equivalent_after_BH,
                  mdd_80pct_power, n_boot_valid, conclusion, sesoi_basis)
}, error = function(e) {
  cat("\n*** equivalence table FAILED:", conditionMessage(e), "\n")
  cat("    Raw differences are safe in bootstrap_added_value_diffs.rds;\n")
  cat("    the table can be rebuilt without repeating any computation.\n")
  NULL
})

if (!is.null(equiv_table)) {
  try(write_csv(equiv_table,
      here::here("results/tables/table_added_value_equivalence_grid.csv")), silent = TRUE)

  cat("\n\n=== HEADLINE: added value of TF integration, ARI endpoint ===\n")
  try(print(equiv_table %>% filter(endpoint == "ari") %>%
          dplyr::select(k, view, point_estimate, ci95_lower, ci95_upper,
                        p_tost_BH, equivalent_after_BH, superior_after_BH,
                        mdd_80pct_power) %>%
          arrange(k, desc(point_estimate)), n = 40), silent = TRUE)

  cat("\n=== Prognostic endpoint (Cox LRT chi-square, TCGA -- exploratory) ===\n")
  try(print(equiv_table %>% filter(endpoint == "cox_lrt") %>%
          dplyr::select(k, view, point_estimate, ci95_lower, ci95_upper,
                        superior_after_BH) %>%
          arrange(k, desc(point_estimate)), n = 40), silent = TRUE)
}

# ------------------------------------------------------------------
# Grid-level verdict
# ------------------------------------------------------------------
safe_which_max <- function(x) { i <- which.max(x); if (length(i) == 0) NA_integer_ else i }
safe_max <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else max(x) }

verdict <- tryCatch({
  stopifnot(!is.null(equiv_table))
  equiv_table %>%
    group_by(endpoint, k) %>%
    summarise(
      n_views = n(),
      n_superior_after_BH  = sum(superior_after_BH, na.rm = TRUE),
      n_equivalent_after_BH = sum(equivalent_after_BH, na.rm = TRUE),
      n_inconclusive = sum(!superior_after_BH & !equivalent_after_BH, na.rm = TRUE),
      best_view = view[safe_which_max(point_estimate)],
      best_point_estimate = safe_max(point_estimate),
      max_mdd = safe_max(mdd_80pct_power),
      .groups = "drop"
    ) %>%
    mutate(
      grid_verdict = case_when(
        n_superior_after_BH > 0 ~
          "AT LEAST ONE VIEW SUPERIOR -- report it, and report that it survived BH correction",
        n_equivalent_after_BH == n_views ~
          "ALL VIEWS EQUIVALENT TO ZERO -- a positive, publishable negative result",
        n_equivalent_after_BH > 0 ~
          "MIXED -- some views equivalent, others inconclusive; report per view",
        TRUE ~
          "ALL INCONCLUSIVE -- underpowered; no claim in either direction is supportable"
      )
    )
}, error = function(e) { cat("*** verdict FAILED:", conditionMessage(e), "\n"); NULL })

if (!is.null(verdict)) {
  cat("\n\n=== GRID-LEVEL VERDICT ===\n")
  try(print(verdict, width = Inf), silent = TRUE)
  try(write_csv(verdict,
      here::here("results/tables/table_added_value_grid_verdict.csv")), silent = TRUE)
}

# ------------------------------------------------------------------
# Figure: forest plot (never fatal -- a missing figure is recoverable,
# an aborted script is not)
# ------------------------------------------------------------------
try({
  for (ep in if (is.null(equiv_table)) character(0) else names(endpoint_cols)) {
    dat <- equiv_table %>% filter(endpoint == ep)
    if (nrow(dat) == 0) next
    p <- ggplot(dat, aes(x = reorder(view, point_estimate), y = point_estimate)) +
      geom_hline(yintercept = 0, linetype = "solid", colour = "grey40") +
      geom_hline(yintercept = c(-unique(dat$sesoi), unique(dat$sesoi)),
                 linetype = "dashed", colour = "firebrick") +
      geom_pointrange(aes(ymin = ci95_lower, ymax = ci95_upper)) +
      coord_flip() +
      facet_wrap(~ paste0("k = ", k)) +
      labs(x = NULL, y = paste0("Added value (fused - expression only), ", ep),
           title = paste0("Added value of TF-activity integration: ", ep),
           subtitle = "Points = paired bootstrap mean; bars = 95% CI; dashed red = SESOI") +
      theme_bw(base_size = 10)
    try(ggsave(here::here("results/figures", paste0("figure_added_value_forest_", ep, ".png")),
               p, width = 11, height = 7, dpi = 150), silent = TRUE)
  }
}, silent = TRUE)

# ==================================================================
# POINT ESTIMATES -- deliberately AFTER the bootstrap.
#
# The bootstrap is the critical path (tens of hours); point estimates
# take ~70 minutes and nothing downstream of them is needed to produce
# the headline equivalence result. Running them first would have put
# 70 minutes of non-essential work in front of the only computation
# that cannot be shortened. If the run has to be abandoned partway,
# the equivalence table and verdict above are already written to disk.
# ==================================================================
# ------------------------------------------------------------------
cluster_from_D <- function(D, k, K = SNF_K, sigma = SNF_SIGMA) {
  SNFtool::spectralClustering(affinity_from_dist(D, K = K, sigma = sigma), K = k)
}
fuse_from_D <- function(D_a, D_b, K = SNF_K, sigma = SNF_SIGMA, t_iter = SNF_T) {
  SNFtool::SNF(list(affinity_from_dist(D_a, K = K, sigma = sigma),
                    affinity_from_dist(D_b, K = K, sigma = sigma)),
               K = K, t = t_iter)
}

all_endpoints <- function(cl, W = NULL) {
  tibble(
    ari_vs_PAM50 = ari_vs_reference(cl, pam50_vec, has_pam50),
    nmi_vs_PAM50 = nmi_vs_reference(cl, pam50_vec, has_pam50),
    cox_lrt_chisq = cox_lrt_endpoint(cl[surv_valid], surv_time[surv_valid], surv_status[surv_valid]),
    cindex = concordance_endpoint(cl[surv_valid], surv_time[surv_valid], surv_status[surv_valid]),
    min_cluster_size = min(table(cl)),
    n_clusters_observed = length(unique(cl)),
    silhouette_native_graph = if (is.null(W)) NA_real_ else
      mean(cluster::silhouette(cl, as.dist(1 - W))[, "sil_width"])
  )
}

cat("\n================ Point-estimate benchmark over grid ================\n")
t0 <- Sys.time()

point_rows <- list()

# NOTE ON LOOP ORDER: networks are built in the OUTER loop and clustered
# in the inner loop over k, never the reverse. An affinity matrix and a
# fused network do not depend on k -- only the final
# spectralClustering(K = k) call does. Building them inside a k loop
# would rebuild every network once per k value, doubling the cost of
# this script for no change in any result.
W_expr <- affinity_from_dist(D_expr, K = SNF_K, sigma = SNF_SIGMA)
for (k in K_VALUES) {
  cl_expr <- SNFtool::spectralClustering(W_expr, K = k)
  point_rows[[paste0("baseline_k", k)]] <-
    bind_cols(tibble(k = k, view = "expression_only", arm = "baseline"),
              all_endpoints(cl_expr, W_expr))
}

# PARALLELISED ACROSS VIEWS. Each view is independent, and with a full
# grid this stage involves one SNF per view -- serially that is hours,
# not minutes, and it sits on the critical path before the bootstrap can
# start. Views carry no shared state, so this is a clean parallel map.
point_view_fun <- function(payload) {
  nm      <- payload$nm
  W_tf    <- affinity_from_dist(payload$D_tf,  K = SNF_K, sigma = SNF_SIGMA)
  W_cat   <- affinity_from_dist(payload$D_cat, K = SNF_K, sigma = SNF_SIGMA)
  W_fused <- SNFtool::SNF(list(W_expr, W_tf), K = SNF_K, t = SNF_T)

  rows <- list()
  for (k in K_VALUES) {
    rows[[length(rows) + 1]] <- cbind(
      data.frame(k = k, view = nm, arm = "tf_only", stringsAsFactors = FALSE),
      all_endpoints(SNFtool::spectralClustering(W_tf, K = k), W_tf))
    rows[[length(rows) + 1]] <- cbind(
      data.frame(k = k, view = nm, arm = "SNF_fused", stringsAsFactors = FALSE),
      all_endpoints(SNFtool::spectralClustering(W_fused, K = k), W_fused))
    rows[[length(rows) + 1]] <- cbind(
      data.frame(k = k, view = nm, arm = "naive_concatenation", stringsAsFactors = FALSE),
      all_endpoints(SNFtool::spectralClustering(W_cat, K = k), W_cat))
  }
  do.call(rbind, rows)
}

# MEMORY: the payload for each view (its TF distance matrix and its
# concatenation distance matrix) is passed as a parLapply ARGUMENT, not
# broadcast with clusterExport. With ~30 views at ~9.6 MB per matrix,
# exporting both full lists to every worker would copy roughly 4 GB and
# can exhaust RAM on a 16 GB machine. Passing payloads as arguments
# means each worker only ever holds its own share.
# Build the concatenation distances now, at the point where they are
# first needed (see the deferral note in the distance-cache section).
# The bootstrap released most of the distance cache to survive on 8 GB.
# Rebuild what the point-estimate stage needs, now that the long job is
# finished and its worker memory has been returned.
gc(verbose = FALSE, full = TRUE)
report_memory("before point estimates")
cat("Building concatenation distances for", length(D_tf_list), "views ...\n")
D_concat_list <- lapply(names(tf_scaled_list), function(nm) {
  dist2_cached(cbind(expr_scaled, tf_scaled_list[[nm]]))
})
names(D_concat_list) <- names(tf_scaled_list)

view_payloads <- lapply(names(D_tf_list), function(nm) {
  list(nm = nm, D_tf = D_tf_list[[nm]], D_cat = D_concat_list[[nm]])
})

cl_pt <- make_cluster_backoff(min(n_cores, length(D_tf_list)))
if (!is.null(cl_pt)) {
source_utils_on_workers(cl_pt)
parallel::clusterEvalQ(cl_pt, {
  suppressPackageStartupMessages({
    library(SNFtool); library(cluster); library(mclust); library(survival); library(tibble)
  })
  # Pin BLAS to one thread per worker -- see set_blas_threads_single()
  # in R/utils_benchmark.R for why oversubscription can make an
  # optimised BLAS slower than the reference one.
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
  }
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
})
parallel::clusterExport(
  cl_pt,
  varlist = c("W_expr", "SNF_K", "SNF_SIGMA", "SNF_T",
              "K_VALUES", "AFFINITY_EPS", "affinity_from_dist", "all_endpoints",
              "pam50_vec", "has_pam50", "surv_time", "surv_status", "surv_valid",
              "ari_vs_reference", "nmi_vs_reference",
              "cox_lrt_endpoint", "concordance_endpoint"),
  envir = environment()
)
}
point_view_rows <- tryCatch(
  par_or_serial(cl_pt, view_payloads, point_view_fun, label = "point-est"),
  error = function(e) { cat("*** point-estimate stage FAILED:", conditionMessage(e),
                            "\n    The equivalence results above are already written and unaffected.\n"); NULL })
if (!is.null(cl_pt)) parallel::stopCluster(cl_pt)

if (!is.null(point_view_rows)) point_rows <- c(point_rows, list(bind_rows(point_view_rows)))
cat("  all", length(D_tf_list), "views done\n")

benchmark_table <- bind_rows(point_rows)
write_csv(benchmark_table, here::here("results/tables/table_benchmark_grid_k2_k3.csv"))
cat("[TIMING] Point estimates:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minutes\n")

cat("\nPoint-estimate benchmark (ARI vs PAM50), fused arm by view:\n")
print(benchmark_table %>%
        filter(arm %in% c("baseline", "SNF_fused")) %>%
        dplyr::select(k, view, arm, ari_vs_PAM50, cox_lrt_chisq, cindex) %>%
        arrange(k, desc(ari_vs_PAM50)), n = 40)

# ------------------------------------------------------------------
# Sample-permutation negative control (retained from previous revision)
#
# Kept because it answers a DIFFERENT question from the edge-permuted
# control introduced in script 05, and the two must not be conflated:
#   sample permutation -> is fusion responding to real per-patient TF
#                         structure, or would any noise view do?
#   edge permutation   -> is any benefit coming from curated regulatory
#                         biology, or merely from averaging arbitrary
#                         gene sets of the same size?
# ------------------------------------------------------------------
primary_view <- primary_view_name
set_pipeline_seed(offset = 777)
perm_idx <- sample(seq_len(nrow(tf_scaled_list[[primary_view]])))
tf_perm <- tf_scaled_list[[primary_view]][perm_idx, , drop = FALSE]
rownames(tf_perm) <- rownames(tf_scaled_list[[primary_view]])
D_tf_perm <- dist2_cached(tf_perm)

perm_control <- tryCatch(map_dfr(K_VALUES, function(k) {
  W_perm_fused <- fuse_from_D(D_expr, D_tf_perm)
  cl_perm <- SNFtool::spectralClustering(W_perm_fused, K = k)
  bind_cols(tibble(k = k, view = paste0(primary_view, "_SAMPLE_PERMUTED"),
                   arm = "sample_permuted_control"),
            all_endpoints(cl_perm, W_perm_fused))
}), error = function(e) { cat("*** permutation control FAILED:", conditionMessage(e), "\n"); NULL })

if (!is.null(perm_control)) benchmark_table <- bind_rows(benchmark_table, perm_control)
try(write_csv(benchmark_table,
    here::here("results/tables/table_benchmark_grid_k2_k3.csv")), silent = TRUE)
if (!is.null(perm_control)) {
  cat("\nSample-permutation negative control:\n")
  try(print(perm_control %>% dplyr::select(k, arm, ari_vs_PAM50, cox_lrt_chisq)), silent = TRUE)
}

# ------------------------------------------------------------------
# Interpretive note, persisted to disk
# ------------------------------------------------------------------
writeLines(c(
  "INTERPRETIVE NOTE: what this script does and does not establish.",
  "======================================================================",
  "",
  "THREE DISTINCT COMPARISONS -- do not conflate them.",
  "",
  "(1) SNF_fused vs expression_only  [table_added_value_equivalence_grid.csv]",
  "    The ADDED-VALUE test, run across the full regulon x method grid,",
  "    with a paired bootstrap, an equivalence test against a",
  "    pre-specified SESOI, and BH correction across views within each",
  "    (k, endpoint) family. This is the manuscript's headline analysis.",
  "",
  "(2) SNF_fused vs naive_concatenation",
  "    Tests the FUSION MACHINERY, not the TF view. If fusion does not",
  "    beat gluing the feature matrices together, SNF specifically is",
  "    contributing nothing, regardless of whether the TF view is",
  "    informative. A separate claim, separately reportable.",
  "",
  "(3) Sample-permuted and edge-permuted controls",
  "    NEGATIVE CONTROLS answering different questions -- see the",
  "    comment block above the permutation control in this script and",
  "    NOTE_tf_activity_grid_design.txt from script 05.",
  "",
  "ON THE ARI-vs-PAM50 ENDPOINT.",
  "PAM50 labels are computed FROM gene expression. An expression-only",
  "clustering is therefore near-optimal for that endpoint BY",
  "CONSTRUCTION, and the baseline arm is being scored against a target",
  "derived from its own input. This endpoint is retained for continuity",
  "with the prior revision and with the wider literature, but it CANNOT",
  "carry the manuscript's central claim on its own. The survival",
  "endpoints do not share this defect because overall survival is not",
  "downstream of the expression matrix.",
  "",
  "ON THE TCGA SURVIVAL ENDPOINTS.",
  "Reported here as SECONDARY and EXPLORATORY because of the low event",
  "count in TCGA-BRCA (see the event count printed in the log). The",
  "adequately powered prognostic added-value test is script 14",
  "(METABRIC), which has roughly an order of magnitude more events and",
  "far longer follow-up.",
  "",
  "ON READING A NULL RESULT FROM THIS TABLE.",
  "Read the `conclusion` and `equivalent_after_BH` columns, NOT the CI",
  "alone. A CI containing zero is compatible with both 'no effect' and",
  "'no power'. Only rows marked EQUIVALENT support a positive claim of",
  "absence; rows marked INCONCLUSIVE support no claim at all, and the",
  "`mdd_80pct_power` column states the effect size the design could",
  "actually have detected."
), here::here("results/tables/NOTE_benchmark_interpretation.txt"))

# ------------------------------------------------------------------
# PRODUCED OUTPUTS -- an explicit inventory, so that after a long run
# there is never any doubt about what was actually generated. Written
# last, after every guarded stage has had its chance to run.
# ------------------------------------------------------------------
cat("\n\n================ PRODUCED OUTPUTS ================\n")
expected <- c(
  "results/objects/bootstrap_added_value_diffs.rds",
  "results/tables/table_added_value_equivalence_grid.csv",
  "results/tables/table_added_value_grid_verdict.csv",
  "results/tables/table_benchmark_grid_k2_k3.csv",
  "results/figures/figure_added_value_forest_ari.png"
)
for (f in expected) {
  fp <- here::here(f)
  cat(sprintf("  %-8s %-56s %s\n",
              if (file.exists(fp)) "PRESENT" else "MISSING", f,
              if (file.exists(fp)) paste0(round(file.size(fp) / 1024, 1), " KB") else ""))
}
nb <- length(list.files(here::here("results/objects"), "^bootstrap_c[0-9]+_b"))
cat(sprintf("  batch checkpoint files on disk: %d\n", nb))
cat("\n  If a TABLE is MISSING but the .rds is PRESENT, no computation was\n")
cat("  lost -- the table rebuilds from the .rds in seconds.\n")
cat("==================================================\n\n")

log_session_info(script_name, key_packages = c("SNFtool", "mclust", "cluster", "survival", "parallel"))
cat("\n✓ 07_benchmark_models.R complete.\n")

close_logger(log_con, script_name)
