# ============================================================
# 19_nested_bootstrap_sensitivity.R   (NEW -- sensitivity analysis)
# Does the primary equivalence verdict survive propagation of
# DISCOVERY-COHORT uncertainty?
# ============================================================
# WHY THIS SCRIPT EXISTS
# ------------------------------------------------------------
# Script 14 establishes the manuscript's primary result: on the
# METABRIC adjusted C-index endpoint, adding a TF-activity view to an
# expression-only SNF clustering is EQUIVALENT TO ZERO against a
# pre-specified SESOI of 0.02, with a minimum detectable difference of
# 0.0012-0.0054 -- four- to seventeenfold below the bound.
#
# That precision was obtained with the TCGA-derived centroids held
# FIXED. Script 14's own header states this and calls a nested
# bootstrap "very large computational cost". It is the single sharpest
# objection available to an examiner, and it is sharp for a specific
# reason: THIS THESIS ARGUES THAT NO NULL MAY BE CLAIMED WITHOUT A
# PRECISION STATEMENT. Leaving the fixed-centroid assumption's effect
# on precision unquantified applies that standard unevenly.
#
# This script quantifies it. Each replicate:
#
#   1. resamples the TCGA discovery cohort with replacement
#   2. RE-DERIVES the expression-only and fused clusterings on that
#      resample (so discovery uncertainty enters the estimate)
#   3. rebuilds centroids from the resampled labels
#   4. projects METABRIC onto them
#   5. also resamples METABRIC patients
#   6. scores the paired added-value difference
#
# Two variance decompositions come from the SAME label sets at no
# extra cost:
#
#   DISCOVERY-ONLY : TCGA resampled, METABRIC fixed.
#                    Isolates the uncertainty script 14 omits.
#   FULL NESTED    : both cohorts resampled.
#                    Total uncertainty. This is the number to report.
#
# HOW TO READ THE OUTPUT -- decide this BEFORE you run it
# ------------------------------------------------------------
# Compare mdd_80pct_power in table_nested_bootstrap_equivalence.csv
# against the pre-specified SESOI of 0.02.
#
#   MDD stays below 0.02  -> the equivalence claim survives EXPLICITLY
#                            rather than by assumption. One paragraph
#                            in Limitations closes the objection.
#
#   MDD rises above 0.02  -> the verdict for these two views becomes
#                            INCONCLUSIVE once discovery uncertainty is
#                            propagated. The point estimates stay near
#                            zero and no view becomes SUPERIOR; what
#                            changes is that the study can no longer
#                            BOUND the effect. Abstract and Conclusion
#                            must then be revised.
#
# Both outcomes are consistent with the thesis. The second requires
# rewriting. Do not start this run unless you are willing to do that.
#
# WHAT IS DELIBERATELY *NOT* PROPAGATED
# ------------------------------------------------------------
# Feature selection is held fixed: the top-2000 and top-5000 variance
# gene sets keep their full-cohort values. Re-selecting features inside
# each replicate would propagate that variance too and would widen
# these intervals further. The MDD produced here is therefore a LOWER
# BOUND on total analytic uncertainty, not an upper one. State this
# rather than letting an examiner find it.
#
# RESAMPLING SCHEME
# ------------------------------------------------------------
# n-out-of-n bootstrap with replacement, matching scripts 07 and 14 so
# the three analyses stay directly comparable. R/utils_benchmark.R's
# make_resample_index() also offers m-out-of-n subsampling, and its own
# comment calls that "the cleaner choice for kNN-graph methods" because
# it cannot create duplicate-sample degeneracy at all. It is NOT used
# here: m-out-of-n subsampling estimates a variance on the wrong scale
# unless rescaled by sqrt(m/n), and silently reporting an unrescaled
# subsample interval alongside script 14's bootstrap interval would
# make the two incomparable -- which is the one thing this script
# exists to avoid. Duplicate-induced degeneracy is instead handled by
# affinity_from_dist()'s epsilon guard, with any replicate that still
# fails recorded and excluded rather than silently dropped.
#
# SCOPE. Only the TWO views for which script 07 computed confidence
# intervals -- dorothea_AC__viper and collectri__viper -- are run here.
# They were fixed before any point estimate was inspected. Running all
# 29 is not affordable and is not necessary: the question is whether
# the precision claim holds, not whether a thirtieth view behaves oddly.
#
# COST -- READ THIS, IT IS THE BINDING CONSTRAINT
# ------------------------------------------------------------
# Dominated by SNFtool::SNF. 07a_timing_probe.R measured one
# view-replicate on this exact cohort and machine at 801.5 core-seconds
# (affinity 3.4 s + SNF 696.5 s + 2 x spectralClustering 14.4 s). SNF is
# 97% of it.
#
# Script 07's own log is the better anchor, because it is measured
# wall-clock rather than modelled: 400 view-replicates on 6 cores took
# 24.78 hours, against a 14.8-hour model prediction. Real throughput was
# therefore ~1300 core-seconds per view-replicate, 1.7x the model. Use
# the measured figure, not the model.
#
# MEASURED ON THIS MACHINE, 2026-08-17. The probe timed two-view
# replicates at 1321 s each, SERIAL. Decomposing against the figures
# above puts SNF at ~634 s per call, consistent with 07a's 696.5 s.
#
# For 98 remaining two-view replicates, after the 8 -> 16 GB upgrade:
#
#     cores  allowance  basis                        wall clock
#       6      1.25x    pre-upgrade projection         ~7.5 h
#       7      1.25x    upgraded, paging removed       ~6.4 h
#       7      1.70x    if contention is not paging    ~8.7 h
#
# BUDGET 6.5-9 HOURS. The spread is honest uncertainty about how much of
# script 07's unexplained 1.7x penalty was memory pressure rather than
# bandwidth contention; the upgrade removes the former but not the
# latter, and which dominated was never measured.
#
# WHY SNF IS SLOW, AND WHY THAT IS NOT BEING FIXED TONIGHT.
# sessionInfo reports "Matrix products: default" -- this pipeline runs
# R's REFERENCE BLAS. SNFtool::SNF is dominated by 1095 x 1095 matrix
# products across 20 cross-diffusion iterations, and an optimised BLAS
# (OpenBLAS/MKL) typically runs those several times faster. The 634 s
# per SNF call is therefore largely a BLAS choice, not an algorithmic
# floor.
#
# It is deliberately left alone. Scripts 07 and 14 produced the
# manuscript's numbers under reference BLAS; computing this sensitivity
# analysis under a different numerical library would mean the check and
# the thing being checked were not computed in the same environment.
# scripts/00b_verify_blas_swap.R exists precisely because that swap has
# to be validated first -- its own header warns that an ILP64 OpenBLAS
# build "does not crash; it silently returns garbage for large
# matrices." That is not a risk to take the night before a defence.
# Worth doing afterwards; not now.
#
# HALVING IT IS LEGITIMATE. Dropping VIEWS_NESTED to the primary view
# alone (dorothea_AC__viper) answers the question the analysis exists to
# answer -- does the precision claim survive discovery uncertainty --
# because that is the view the manuscript's primary analysis uses. The
# second view is corroboration, not the test. If time is short, run one
# view and say so.
#
# Distance matrices are CACHED and indexed, exactly as in
# utils_benchmark.R, so the O(n^2 p) step is paid once rather than per
# replicate. RUN THE PROBE FIRST -- it reports a projected wall-clock
# total from YOUR machine before you commit, and its replicates are
# cached, so probing costs nothing the full run would not have paid.
#
# IF THIS SCRIPT STOPS EARLY
# ------------------------------------------------------------
# init_logger() redirects console output with sink(). If the script
# halts before close_logger() runs, your R console will appear silent.
# Restore it with:
#
#     sink(type = "message"); sink()
#
# Every deliberate halt below calls fail_clean(), which restores the
# sinks first. Only an unanticipated error can leave them open.
# ============================================================

# MEMORY: org.Hs.eg.db and survival are deliberately NOT loaded here.
# They are needed only in stage B, which runs after the cluster has been
# stopped. org.Hs.eg.db is a large SQLite-backed annotation package;
# holding it in the master for the whole of stage A costs hundreds of MB
# at exactly the moment six workers are competing for the same 8 GB.
# R/utils.R and R/utils_benchmark.R require only `here` and `tidyverse`,
# and map_to_symbol_dedup() reaches AnnotationDbi through `::`, so
# deferring is safe. They are loaded at the top of stage B.
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SNFtool)
})
source(here::here("R", "utils.R"))
source(here::here("R", "utils_benchmark.R"))

script_name <- "19_nested_bootstrap_sensitivity"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

# Restore console sinks before halting, so a deliberate stop does not
# leave the user's session mute.
fail_clean <- function(...) {
  msg <- paste0(...)
  try(close_logger(log_con, script_name), silent = TRUE)
  stop(msg, call. = FALSE)
}

# Memory is the binding constraint on this machine (8 GB total), and
# script 07 lost two runs to it. This reports only what the MASTER R
# process holds -- worker memory is on top of it, roughly 320 MB each.
# Watch these lines: if the master climbs unexpectedly during stage A,
# something is being retained that should have been rm()'d.
report_memory <- function(tag) {
  g <- gc(verbose = FALSE)
  cat(sprintf("[memory %s] master R is holding ~%.0f MB\n", tag, sum(g[, 2])))
  utils::flush.console()
}

# ------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------
# TRUE  -> run N_PROBE replicates, print a projected total, stop.
#          ALWAYS do this first. Probe replicates are cached and
#          reused by the full run, so nothing is wasted.
PROBE_ONLY <- FALSE   # probe already passed 2026-08-17; set TRUE to re-probe
N_PROBE    <- 2

N_BOOT_NESTED <- 100      # 100 is defensible; below 50 TOST refuses

# CORE COUNT -- REVISED 2026-08-17 AFTER THE RAM UPGRADE (8 GB -> 16 GB)
#
# The old ceiling of six was never about CPU. Script 07's header records
# why it existed:
#
#   "Each PSOCK worker costs roughly 320 MB [...] and the master holds
#    a further ~800 MB. Seven workers is ~3.0 GB of R on top of the
#    3-4 GB Windows and RStudio already use -- which is why a worker
#    was killed 7.5 hours into the previous attempt ('error reading
#    from connection'). [...] this run has already been lost twice to
#    memory pressure."
#
# That constraint is gone. At 16 GB, seven workers is ~2.3 GB plus a
# master of roughly 1.0-1.2 GB (lower than before, since org.Hs.eg.db
# and survival are no longer loaded during stage A), against 3-4 GB for
# Windows and RStudio: about 7 GB of 16, with real headroom rather than
# none. Seven is now the right number, and the cap at
# detectCores() - 1 keeps one core free for the master and the OS.
#
# A SECOND, LESS OBVIOUS GAIN. Script 07 ran 1.7x slower than its
# core-second model predicted, and that gap was never fully explained.
# Some of it is dispatch and memory-bandwidth contention, but some was
# very likely paging: at 8 GB the machine was operating with almost no
# free memory, and Windows pages rather than failing. With 16 GB that
# component should disappear, so the realistic allowance moves back
# toward the 1.25x in the projection rather than script 07's 1.7x.
N_CORES_19    <- 7

# ANYTHING ELSE RUNNING ON THIS MACHINE?
# At 16 GB this is no longer a binding concern: seven workers plus the
# master come to roughly 3.5 GB, leaving several GB free even with a
# browser open. The [memory] lines below still report the master's
# footprint -- watch them out of habit, not anxiety.

K_VALUES  <- c(2, 3)
SNF_K     <- 20
SNF_SIGMA <- 0.5
SNF_T     <- 20

SESOI_CINDEX     <- 0.02   # UNCHANGED from script 14. Do not edit.
MIN_COVERAGE_PCT <- 70
MIN_CLUSTER_N    <- 5      # centroid needs at least this many members

# The two views script 07 bootstrapped, fixed before any point estimate
# was examined. Do not EXTEND this list after seeing results.
#
# Dropping to the first view alone roughly halves the runtime and is
# legitimate -- dorothea_AC__viper is the manuscript's primary view, so
# it carries the test; collectri__viper is corroboration. Reducing scope
# for time is a reporting decision, not a forking path, provided you say
# which views were run. Choose BEFORE the run, not after seeing a result
# you dislike.
VIEWS_NESTED <- c("dorothea_AC__viper", "collectri__viper")
# VIEWS_NESTED <- c("dorothea_AC__viper")   # <- uncomment to halve runtime

# ------------------------------------------------------------------
# STAGE 0. Inputs (deliberately WITHOUT METABRIC -- stage A must not
# hold the ~300 MB expression matrix while five workers are running)
# ------------------------------------------------------------------
vst_top2000     <- readRDS(here::here("data/processed/vst_top2000_genes.rds"))
vst_top5000     <- readRDS(here::here("data/processed/vst_top5000_genes.rds"))
tf_grid         <- readRDS(here::here("data/processed/tf_activity_grid.rds"))
sample_metadata <- readRDS(here::here("data/processed/sample_metadata_matched.rds"))

missing_views <- setdiff(VIEWS_NESTED, names(tf_grid))
if (length(missing_views)) {
  fail_clean("Views not present in tf_activity_grid.rds: ",
             paste(missing_views, collapse = ", "))
}
tf_grid <- tf_grid[VIEWS_NESTED]

common_samples <- Reduce(intersect, c(
  list(colnames(vst_top2000), colnames(vst_top5000), rownames(sample_metadata)),
  lapply(tf_grid, colnames)
))
vst_top2000 <- vst_top2000[, common_samples, drop = FALSE]
vst_top5000 <- vst_top5000[, common_samples, drop = FALSE]
tf_grid     <- lapply(tf_grid, function(m) m[, common_samples, drop = FALSE])
n_tcga <- length(common_samples)
cat("TCGA discovery samples:", n_tcga, "\n")
cat("Views:", paste(VIEWS_NESTED, collapse = ", "), "\n")

# ------------------------------------------------------------------
# STAGE A. Cached distances, verified against SNFtool::dist2
#
# d(i,j) is a property of the PAIR and does not change under
# resampling, so each replicate indexes D[idx, idx] rather than
# recomputing. verify_distance_cache() FAILS HARD on any convention
# mismatch -- a fast wrong answer here would silently invalidate the
# whole sensitivity analysis.
#
# NOTE: verify_distance_cache() calls set_pipeline_seed(offset = 9001)
# internally, so the RNG state after this block is not the state set at
# the top of the script. That is harmless because every replicate seeds
# itself deterministically from its own index.
# ------------------------------------------------------------------
scale_view <- function(mat) { s <- scale(t(mat)); s[is.na(s)] <- 0; s }

expr_scaled <- scale_view(vst_top2000)
D_expr <- dist2_cached(expr_scaled)
verify_distance_cache(expr_scaled, D_expr, label = "expression (nested)")

D_tf <- list()
for (nm in VIEWS_NESTED) {
  s_tf <- scale_view(tf_grid[[nm]])
  D_tf[[nm]] <- dist2_cached(s_tf)
  verify_distance_cache(s_tf, D_tf[[nm]], label = paste0("TF: ", nm))
  rm(s_tf)
}
rm(expr_scaled, tf_grid, vst_top2000, sample_metadata)
gc(verbose = FALSE, full = TRUE)
cat("Distance caches built and verified.\n")
report_memory("after distance caches")

# ------------------------------------------------------------------
# PREFLIGHT -- verify every STAGE B prerequisite BEFORE stage A starts.
#
# Stage A is ~7.5 hours. Stage B needs files and columns that stage A
# never touches, so without this block a missing input would surface
# only after the whole overnight run had completed. Discovering at 8 a.m.
# that a five-second export was never made is not an acceptable failure
# mode. Everything below is checked in about ten seconds.
# ------------------------------------------------------------------
cat("\n--- preflight: stage B prerequisites ---\n")

pf_files <- c(
  "METABRIC expression"   = here::here("data/external/metabric_expression.txt"),
  "METABRIC clinical"     = here::here("data/external/metabric_clinical.txt"),
  "All-gene VST matrix"   = here::here("data/processed/vst_matrix_all_genes.rds"),
  "Prepared clinical frame (from scripts/14b_export_clin.R)" =
    here::here("results/objects/metabric_clin_prepared.rds")
)
pf_missing <- pf_files[!file.exists(pf_files)]
if (length(pf_missing)) {
  fail_clean("Stage B inputs missing -- fix these BEFORE starting the run:\n",
             paste(sprintf("  [%s]\n    %s", names(pf_missing), pf_missing),
                   collapse = "\n"),
             "\n\nIf the prepared clinical frame is the missing one, run:\n",
             "  source(here::here('scripts/14b_export_clin.R'))   # ~5 seconds")
}
for (nm in names(pf_files)) cat("  OK  ", nm, "\n")

# The clinical frame must carry every column the Cox models use.
pf_clin <- readRDS(pf_files[["Prepared clinical frame (from scripts/14b_export_clin.R)"]])
pf_need <- c("sample_id", "os_years", "os_status", "age_covariate", "er_covariate")
if (!all(pf_need %in% names(pf_clin))) {
  fail_clean("Prepared clinical frame is missing columns: ",
             paste(setdiff(pf_need, names(pf_clin)), collapse = ", "),
             "\nRe-run scripts/14b_export_clin.R.")
}
cat("  OK   clinical frame carries all five model variables (",
    nrow(pf_clin), "rows )\n")

# tcga_z is built from vst_top5000 in stage B and must line up column-wise
# with common_samples, since centroids are indexed by the resample index.
# Checking the precondition here avoids discovering it after stage A.
if (!identical(colnames(vst_top5000), common_samples)) {
  fail_clean("colnames(vst_top5000) does not match common_samples. ",
             "Centroid construction in stage B indexes by position and ",
             "would be silently wrong.")
}
cat("  OK   vst_top5000 columns align with common_samples\n")
rm(pf_clin, pf_files, pf_missing, pf_need); gc(verbose = FALSE, full = TRUE)
cat("Preflight passed -- stage B will not fail for a missing input.\n")

# ------------------------------------------------------------------
# One replicate of the DISCOVERY step.
#
# ERROR CONTAINMENT IS LOAD-BEARING HERE. affinity_from_dist() and
# spectralClustering() both signal errors rather than returning a
# sentinel; inside parLapply a single signalling replicate aborts the
# ENTIRE chunk and the script with it, discarding hours of completed
# SNF. The whole body is therefore wrapped so that a failed replicate
# returns ok = FALSE and is counted and excluded, never silently
# dropped and never fatal.
#
# ON LABEL SWITCHING. Spectral cluster IDs are arbitrary and are NOT
# comparable across replicates. They do not need to be: centroids are
# built from within-replicate labels, METABRIC is classified by nearest
# centroid within the same replicate, and the Cox model treats the
# assignment as an unordered factor. No cross-replicate label matching
# is required, and adding one would be a bug.
# ------------------------------------------------------------------
nested_replicate <- function(b) {
  tryCatch({
    set.seed(PIPELINE_SEED + 19000L + b)
    idx <- sample.int(n_tcga, n_tcga, replace = TRUE)

    We <- affinity_from_dist(D_expr[idx, idx, drop = FALSE],
                             K = SNF_K, sigma = SNF_SIGMA)

    labels <- list()
    for (k in K_VALUES) {
      labels[[paste0("expression_only|k", k)]] <-
        SNFtool::spectralClustering(We, K = k)
    }

    for (nm in names(D_tf)) {
      Wt <- affinity_from_dist(D_tf[[nm]][idx, idx, drop = FALSE],
                               K = SNF_K, sigma = SNF_SIGMA)
      Wf <- SNFtool::SNF(list(We, Wt), K = SNF_K, t = SNF_T)
      if (!all(is.finite(Wf))) {
        stop("non-finite fused network for view ", nm)
      }
      for (k in K_VALUES) {
        labels[[paste0(nm, "|k", k)]] <- SNFtool::spectralClustering(Wf, K = k)
      }
      rm(Wt, Wf)
    }

    list(b = b, ok = TRUE, idx = idx, labels = labels, reason = NA_character_)
  },
  error = function(e) {
    list(b = b, ok = FALSE, idx = NULL, labels = NULL,
         reason = conditionMessage(e))
  })
}

# Per-replicate checkpointing. Finer than batch-level: the probe's
# replicates are reused by the full run, and an interrupted run resumes
# at the exact replicate it reached.
rep_path <- function(b) {
  here::here("results/objects", sprintf("nested_rep_%04d.rds", b))
}

# ------------------------------------------------------------------
# TIMING PROBE
# ------------------------------------------------------------------
cat("\n--- timing probe:", N_PROBE, "replicate(s), serial ---\n")
utils::flush.console()

n_cached <- 0L
t_probe <- Sys.time()
probe <- lapply(seq_len(N_PROBE), function(b) {
  if (file.exists(rep_path(b))) {
    cat("  replicate", b, "already cached\n")
    n_cached <<- n_cached + 1L
    return(readRDS(rep_path(b)))
  }
  r <- nested_replicate(b)
  saveRDS(r, rep_path(b))
  r
})
n_timed <- N_PROBE - n_cached
probe_secs <- if (n_timed > 0) {
  as.numeric(difftime(Sys.time(), t_probe, units = "secs")) / n_timed
} else NA_real_

failed <- !vapply(probe, function(x) isTRUE(x$ok), logical(1))
if (any(failed)) {
  fail_clean("PROBE FAILED: ",
             paste(vapply(probe[failed], function(x) x$reason, character(1)),
                   collapse = "; "))
}
if (is.na(probe_secs)) {
  cat("Probe PASS (all probe replicates were already cached, so no timing\n",
      "was measured this run -- the projection below is skipped).\n", sep = "")
} else {
  cat(sprintf("Probe PASS. %.0f s per replicate (serial, %d timed).\n",
              probe_secs, n_timed))
}
cat(sprintf("  k=2 sizes, expression arm : %s\n",
            paste(table(probe[[1]]$labels[["expression_only|k2"]]), collapse = "/")))
for (nm in VIEWS_NESTED) {
  cat(sprintf("  k=2 sizes, %-22s: %s\n", nm,
              paste(table(probe[[1]]$labels[[paste0(nm, "|k2")]]), collapse = "/")))
}

n_cores_19 <- max(1, min(N_CORES_19, parallel::detectCores() - 1))
# 1.25x allowance for parallel dispatch overhead and memory-bandwidth
# contention: five workers do not deliver a clean fivefold speed-up on
# matrices this size.
if (!is.na(probe_secs)) {
  n_remaining <- sum(!file.exists(vapply(seq_len(N_BOOT_NESTED), rep_path, character(1))))
  proj_hours <- probe_secs * n_remaining * 1.25 / n_cores_19 / 3600
  cat(sprintf(
    "\nPROJECTED STAGE-A WALL CLOCK: %.1f hours for the %d remaining replicates on %d cores.\n",
    proj_hours, n_remaining, n_cores_19))
  cat("  (includes a 1.25x allowance for dispatch and memory-bandwidth contention;\n")
  cat("   cached replicates are excluded)\n")
  cat("Stage B (centroids + projection + Cox) adds roughly 10-20 minutes.\n")
}

if (PROBE_ONLY) {
  cat("\nPROBE_ONLY is TRUE -- stopping here, by design.\n")
  cat("Affordable?  set PROBE_ONLY <- FALSE and re-run.\n")
  cat("Too slow?    lower N_BOOT_NESTED (not below 60) and re-run this probe.\n")
  fail_clean("Probe complete -- this halt is intentional, not an error.")
}

# ------------------------------------------------------------------
# STAGE A (full). Resume-safe: replicates already on disk are skipped.
# ------------------------------------------------------------------
all_b <- seq_len(N_BOOT_NESTED)
todo  <- all_b[!file.exists(vapply(all_b, rep_path, character(1)))]
cat(sprintf("\n--- stage A: %d of %d replicates still to run ---\n",
            length(todo), N_BOOT_NESTED))
cat("\n  !! THIS RUN IS UNATTENDED AND LONG. Before you walk away:\n")
cat("     - disable Windows sleep/hibernate (Settings > Power > Screen and sleep\n")
cat("       > 'When plugged in, put my device to sleep after' = Never)\n")
cat("     - keep the laptop on mains power\n")
cat("     - at 16 GB with 7 workers there is ample memory headroom\n")
cat("     A sleeping machine suspends the workers and the run does not finish.\n")
cat("     Completed replicates are on disk, so re-sourcing resumes -- but you\n")
cat("     lose the night.\n\n")
utils::flush.console()

if (length(todo)) {
  cl19 <- tryCatch(parallel::makeCluster(n_cores_19), error = function(e) NULL)
  if (!is.null(cl19)) {
    parallel::clusterEvalQ(cl19, {
      suppressPackageStartupMessages(library(SNFtool))
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
        try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
      }
      Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
                 MKL_NUM_THREADS = "1")
    })
    # affinity_from_dist() closes over AFFINITY_EPS in the global
    # environment, so AFFINITY_EPS must be exported alongside it.
    parallel::clusterExport(
      cl19,
      varlist = c("D_expr", "D_tf", "n_tcga", "SNF_K", "SNF_SIGMA", "SNF_T",
                  "K_VALUES", "PIPELINE_SEED", "AFFINITY_EPS",
                  "affinity_from_dist", "nested_replicate"),
      envir = environment()
    )
  } else {
    cat("  [warn] could not start a cluster; running serially.\n")
  }

  # CHUNK SIZE = one round per worker, not two. Results are written only
  # after a whole chunk returns, so the chunk size IS the crash-exposure
  # window. At ~22 minutes per replicate, chunks of 2*n_cores would put
  # 44 minutes of compute at risk; one round halves that to ~22 minutes.
  # The extra cluster round-trips cost seconds against chunks that take
  # twenty minutes.
  chunks <- split(todo, ceiling(seq_along(todo) / n_cores_19))
  for (ci in seq_along(chunks)) {
    ch <- chunks[[ci]]
    t_c <- Sys.time()
    res <- if (!is.null(cl19)) {
      tryCatch(parallel::parLapply(cl19, ch, nested_replicate),
               error = function(e) { cat("  [warn] chunk failed on cluster (",
                 conditionMessage(e), ") -- retrying serially\n", sep = "")
                 lapply(ch, nested_replicate) })
    } else {
      lapply(ch, nested_replicate)
    }
    for (r in res) saveRDS(r, rep_path(r$b))
    n_bad <- sum(!vapply(res, function(x) isTRUE(x$ok), logical(1)))
    cat(sprintf("  chunk %d/%d (reps %d-%d) done in %.1f min | %d failed\n",
                ci, length(chunks), min(ch), max(ch),
                as.numeric(difftime(Sys.time(), t_c, units = "mins")), n_bad))
    if (ci %% 3 == 1) report_memory(sprintf("after chunk %d", ci))
    utils::flush.console()
  }
  if (!is.null(cl19)) try(parallel::stopCluster(cl19), silent = TRUE)
}

nested_labels <- lapply(all_b, function(b) {
  if (file.exists(rep_path(b))) readRDS(rep_path(b)) else NULL
})
nested_labels <- Filter(Negate(is.null), nested_labels)
ok_flag <- vapply(nested_labels, function(x) isTRUE(x$ok), logical(1))
if (any(!ok_flag)) {
  cat("Failed replicates and reasons:\n")
  for (r in nested_labels[!ok_flag]) cat("  b =", r$b, ":", r$reason, "\n")
}
cat(sprintf("Stage A complete: %d/%d replicates valid.\n",
            sum(ok_flag), length(ok_flag)))
if (sum(ok_flag) < 50) {
  fail_clean("Fewer than 50 valid replicates; tost_equivalence_bootstrap() ",
             "will refuse to run and no equivalence claim would be meaningful.")
}
nested_labels <- nested_labels[ok_flag]

rm(D_expr, D_tf); gc(verbose = FALSE, full = TRUE)

# ------------------------------------------------------------------
# STAGE B. Centroids, METABRIC projection, scoring.
#
# Machinery is IDENTICAL to script 14 so the two analyses differ only
# in what is resampled. Any divergence here would confound the
# comparison this script exists to make.
# ------------------------------------------------------------------
cat("\n--- stage B: METABRIC projection and scoring ---\n")
report_memory("stage B start -- workers released")

# Loaded HERE, not at the top of the script: see the memory note beside
# the library() block. The cluster is stopped by this point, so the
# master is free to take the space.
suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(survival)
})

vst_matrix_all_genes <- readRDS(here::here("data/processed/vst_matrix_all_genes.rds"))
mapped_5000 <- map_to_symbol_dedup(vst_top5000, variance_reference = vst_matrix_all_genes,
                                   org_db = org.Hs.eg.db)
tcga_z <- t(scale(t(mapped_5000$matrix)))
rm(vst_matrix_all_genes, mapped_5000); gc(verbose = FALSE, full = TRUE)
cat("TCGA centroid feature space:", nrow(tcga_z), "gene symbols\n")
stopifnot(identical(colnames(tcga_z), common_samples))

metabric_expr_path <- here::here("data/external/metabric_expression.txt")
metabric_clin_path <- here::here("data/external/metabric_clinical.txt")
if (!file.exists(metabric_expr_path)) fail_clean("METABRIC expression file not found.")

metabric_raw <- read.delim(metabric_expr_path, check.names = FALSE)
gene_col <- intersect(c("Hugo_Symbol", "GENE_SYMBOL", "gene_symbol"),
                      colnames(metabric_raw))[1]
if (is.na(gene_col)) fail_clean("No gene-symbol column in METABRIC expression file.")
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

common_genes <- intersect(rownames(tcga_z), rownames(metabric_expr))
coverage_pct <- 100 * length(common_genes) / nrow(tcga_z)
cat(sprintf("Gene coverage: %d/%d (%.1f%%)\n",
            length(common_genes), nrow(tcga_z), coverage_pct))
if (coverage_pct < MIN_COVERAGE_PCT) {
  fail_clean("Gene coverage below ", MIN_COVERAGE_PCT, "%.")
}

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
rm(metabric_expr, metabric_common); gc(verbose = FALSE, full = TRUE)
cat("Final common gene set:", length(common_genes), "\n")

# --- METABRIC clinical: reuse script 14's prepared frame ---
clin_cache <- here::here("results/objects/metabric_clin_prepared.rds")
if (!file.exists(clin_cache)) {
  fail_clean(
    "Expected ", clin_cache, ".\n",
    "Script 14 builds the METABRIC clinical frame (`clin`) inline. Add this line\n",
    "to scripts/14_prognostic_added_value.R immediately after the\n",
    "  cat(\"METABRIC OS status resolved:\\n\"); print(table(clin$os_status, ...))\n",
    "block, re-run script 14 up to that point, then re-run this script:\n\n",
    "  saveRDS(clin, here::here('results/objects/metabric_clin_prepared.rds'))\n\n",
    "Re-deriving the clinical recoding here instead would risk the two analyses\n",
    "silently diverging on OS-status or ER coding -- the ER column contains the\n",
    "cBioPortal misspelling 'Positve', which script 14 handles explicitly. A\n",
    "divergence there is exactly the confound this sensitivity analysis must\n",
    "not introduce."
  )
}
clin <- readRDS(clin_cache)
required_clin <- c("sample_id", "os_years", "os_status", "age_covariate", "er_covariate")
if (!all(required_clin %in% names(clin))) {
  fail_clean("Cached clinical frame is missing: ",
             paste(setdiff(required_clin, names(clin)), collapse = ", "))
}

base_df <- tibble(sample_id = colnames(metabric_z)) %>%
  left_join(clin, by = "sample_id") %>%
  filter(is.finite(os_years), os_years > 0, !is.na(os_status),
         !is.na(age_covariate), !is.na(er_covariate))
cat("METABRIC analysis set:", nrow(base_df), "patients,",
    sum(base_df$os_status == 1), "events\n")
cat("  (script 14 reported 1937 patients / 1125 events -- these must match)\n")

# `sub` is the resampled TCGA expression block. It depends only on idx,
# which is constant across the six arms within a replicate, so it is
# built ONCE per replicate by the caller rather than six times here.
# Each build copies ~4545 x 1095 doubles (~40 MB); doing it per arm cost
# five redundant 40 MB copies per replicate.
build_centroids_sub <- function(labels_vec, sub) {
  grp <- split(seq_len(ncol(sub)), paste0("SNF_C", labels_vec))
  # A cluster too small to define a stable centroid is dropped rather
  # than producing a degenerate centroid and a spurious METABRIC
  # assignment. Returning NULL propagates as an honest missing value.
  grp <- grp[vapply(grp, length, integer(1)) >= MIN_CLUSTER_N]
  if (length(grp) < 2) return(NULL)
  sapply(grp, function(s) rowMeans(sub[, s, drop = FALSE]))
}

# `mz_sub` is metabric_z already restricted to the common gene set. That
# subset is constant, so it is taken ONCE outside the replicate loop
# instead of copying ~3884 x 1980 doubles (~60 MB) on every one of the
# six calls per replicate.
classify_to_centroids <- function(mz_sub, centroids, genes) {
  cmat <- centroids[genes, , drop = FALSE]
  cor_matrix <- suppressWarnings(cor(mz_sub, cmat, method = "pearson"))
  # A zero-variance centroid yields NA correlations; max.col would then
  # pick column 1 for every patient and manufacture a partition.
  if (anyNA(cor_matrix)) return(NULL)
  best_idx <- max.col(cor_matrix, ties.method = "first")
  setNames(colnames(cmat)[best_idx], rownames(cor_matrix))
}

score_arm <- function(cluster_vec, df, rows) {
  d <- df[rows, , drop = FALSE]
  d$cl <- factor(cluster_vec[d$sample_id])
  if (length(unique(d$cl)) < 2) {
    return(c(cindex_adj = NA_real_, cindex_cl = NA_real_))
  }
  fit_full <- try(coxph(Surv(os_years, os_status) ~ cl + age_covariate + er_covariate,
                        data = d), silent = TRUE)
  fit_cl   <- try(coxph(Surv(os_years, os_status) ~ cl, data = d), silent = TRUE)
  if (inherits(fit_full, "try-error")) {
    return(c(cindex_adj = NA_real_, cindex_cl = NA_real_))
  }
  c(cindex_adj = unname(summary(fit_full)$concordance[1]),
    cindex_cl  = if (inherits(fit_cl, "try-error")) NA_real_
                 else unname(summary(fit_cl)$concordance[1]))
}

all_rows <- seq_len(nrow(base_df))

# Hoisted out of the replicate loop: constant for the whole run.
metabric_z_sub <- metabric_z[common_genes, , drop = FALSE]

nested_rows <- list()
n_skipped <- 0L
for (i in seq_along(nested_labels)) {
  rep_i <- nested_labels[[i]]
  b <- rep_i$b; idx <- rep_i$idx

  # Constant across the six arms of this replicate -- built once here.
  sub_tcga <- tcga_z[, common_samples[idx], drop = FALSE]

  # METABRIC resample for this replicate. The SAME rows are used for
  # every arm within the replicate, so the difference stays a paired
  # contrast under both schemes.
  set.seed(PIPELINE_SEED + 19500L + b)
  rows_mb <- sample(all_rows, length(all_rows), replace = TRUE)

  arms <- names(rep_i$labels)
  assign_list <- vector("list", length(arms))
  names(assign_list) <- arms
  for (arm in arms) {
    cents <- build_centroids_sub(rep_i$labels[[arm]], sub_tcga)
    # NOTE: `assign_list[[arm]] <- NULL` would DELETE the element rather
    # than store a NULL, so the assignment is made via list() to keep
    # the slot present and the lookup below unambiguous.
    assign_list[arm] <- list(
      if (is.null(cents)) NULL else
        classify_to_centroids(metabric_z_sub, cents, common_genes)
    )
  }

  for (k in K_VALUES) {
    base_arm <- paste0("expression_only|k", k)
    if (is.null(assign_list[[base_arm]])) { n_skipped <- n_skipped + 1L; next }
    base_fix  <- score_arm(assign_list[[base_arm]], base_df, all_rows)
    base_full <- score_arm(assign_list[[base_arm]], base_df, rows_mb)

    for (nm in VIEWS_NESTED) {
      a <- assign_list[[paste0(nm, "|k", k)]]
      if (is.null(a)) { n_skipped <- n_skipped + 1L; next }
      s_fix  <- score_arm(a, base_df, all_rows)
      s_full <- score_arm(a, base_df, rows_mb)
      nested_rows[[length(nested_rows) + 1]] <- data.frame(
        b = b, k = k, view = nm,
        d_cindex_adj_discovery_only = s_fix[["cindex_adj"]]  - base_fix[["cindex_adj"]],
        d_cindex_cl_discovery_only  = s_fix[["cindex_cl"]]   - base_fix[["cindex_cl"]],
        d_cindex_adj_full_nested    = s_full[["cindex_adj"]] - base_full[["cindex_adj"]],
        d_cindex_cl_full_nested     = s_full[["cindex_cl"]]  - base_full[["cindex_cl"]]
      )
    }
  }
  if (i %% 10 == 0) {
    cat(sprintf("  scored %d/%d replicates\n", i, length(nested_labels)))
    utils::flush.console()
  }
}

if (!length(nested_rows)) fail_clean("No replicate produced a scoreable contrast.")
boot_nested <- bind_rows(nested_rows) %>% as_tibble()
saveRDS(boot_nested, here::here("results/objects/bootstrap_nested_discovery.rds"))
cat("Scored replicates:", length(unique(boot_nested$b)),
    "| arm-k combinations skipped for degenerate centroids:", n_skipped, "\n")

# ------------------------------------------------------------------
# STAGE C. Equivalence testing and the comparison that matters
#
# BH is applied WITHIN each (k, endpoint) family, matching script 14
# exactly. Note that each family here contains only 2 views rather than
# 29, so the correction is nearly inert; it is retained for procedural
# consistency, not because it does much work.
# ------------------------------------------------------------------
ep_map <- c(
  cindex_adjusted_full_nested     = "d_cindex_adj_full_nested",
  cindex_adjusted_discovery_only  = "d_cindex_adj_discovery_only",
  cindex_cluster_only_full_nested = "d_cindex_cl_full_nested"
)

equiv_nested <- tryCatch({
  map_dfr(names(ep_map), function(ep) {
    col <- ep_map[[ep]]
    map_dfr(K_VALUES, function(kk) {
      map_dfr(VIEWS_NESTED, function(vv) {
        d <- boot_nested[[col]][boot_nested$k == kk & boot_nested$view == vv]
        d <- d[is.finite(d)]
        if (length(d) < 50) {
          cat("  [skip]", ep, vv, paste0("k", kk), "-- only", length(d),
              "finite replicates\n")
          return(tibble())
        }
        tost_equivalence_bootstrap(d, sesoi = SESOI_CINDEX,
                                   label = paste(ep, vv, paste0("k", kk))) %>%
          mutate(endpoint = ep, k = kk, view = vv,
                 mdd_80pct_power = minimum_detectable_difference(d),
                 boot_sd = sd(d))
      })
    })
  }) %>%
    group_by(k, endpoint) %>%
    group_modify(~ adjust_grid_multiplicity(.x)) %>%
    ungroup()
}, error = function(e) {
  cat("*** equivalence table FAILED:", conditionMessage(e),
      "\n    Raw bootstrap is safe in results/objects/ -- rebuildable in seconds.\n")
  NULL
})

if (is.null(equiv_nested) || !nrow(equiv_nested)) {
  fail_clean("Equivalence table is empty; inspect bootstrap_nested_discovery.rds.")
}
write_csv(equiv_nested, here::here("results/tables/table_nested_bootstrap_equivalence.csv"))

# --- side-by-side against script 14 ---
fixed_path <- here::here("results/tables/table_metabric_prognostic_added_value.csv")
if (file.exists(fixed_path)) {
  fixed <- read_csv(fixed_path, show_col_types = FALSE)
  fixed_sub <- fixed %>%
    filter(endpoint == "cindex_adjusted", view %in% VIEWS_NESTED) %>%
    transmute(view, k, scheme = "fixed centroids (script 14)",
              point_estimate, mdd_80pct_power, conclusion)
  nested_sub <- equiv_nested %>%
    filter(endpoint == "cindex_adjusted_full_nested") %>%
    transmute(view, k, scheme = "full nested (script 19)",
              point_estimate, mdd_80pct_power, conclusion)
  comparison <- bind_rows(fixed_sub, nested_sub) %>% arrange(view, k, scheme)
  write_csv(comparison, here::here("results/tables/table_nested_vs_fixed_precision.csv"))
  cat("\n=== PRECISION UNDER THE TWO BOOTSTRAP SCHEMES ===\n")
  print(comparison, n = 40, width = Inf)
} else {
  cat("\n[warn] script 14's added-value table not found; side-by-side skipped.\n")
}

cat("\n=== NESTED BOOTSTRAP, primary endpoint (full nested) ===\n")
# dplyr::select(), NOT bare select(). library(org.Hs.eg.db) above pulls in
# AnnotationDbi, whose select() is an S4 generic that MASKS dplyr's and has
# no method for a tbl_df. The 2026-08-17 run crashed here after 6.5 hours of
# completed stage A/B with:
#   Error: unable to find an inherited method for function 'select'
#   for signature 'x = "tbl_df"'
# Script 14 avoids the same trap by namespacing its select() calls; do not
# un-namespace this one. Every other dplyr verb used below this point
# (filter, transmute, bind_rows, arrange, mutate) was exercised under the
# same masking in that run and is unaffected -- select is the only casualty.
print(equiv_nested %>% dplyr::filter(endpoint == "cindex_adjusted_full_nested") %>%
        dplyr::select(view, k, point_estimate, ci90_lower, ci90_upper,
                      boot_sd, mdd_80pct_power, equivalent_after_BH, conclusion),
      n = 20, width = Inf)

primary <- equiv_nested %>% filter(endpoint == "cindex_adjusted_full_nested")
n_tested <- nrow(primary)
if (n_tested == 0) {
  verdict <- "NO PRIMARY-ENDPOINT ROWS PRODUCED -- inspect the raw bootstrap."
} else {
  max_mdd  <- suppressWarnings(max(primary$mdd_80pct_power, na.rm = TRUE))
  n_equiv  <- sum(primary$equivalent_after_BH, na.rm = TRUE)
  n_sup    <- sum(primary$superior_after_BH, na.rm = TRUE)
  verdict <- if (is.finite(max_mdd) && max_mdd < SESOI_CINDEX && n_equiv == n_tested) {
    sprintf(paste0(
      "EQUIVALENCE SURVIVES. Under the full nested bootstrap the largest MDD is ",
      "%.4f, still below the pre-specified SESOI of %.2f, and %d/%d view-by-k ",
      "combinations remain equivalent after BH. The primary claim no longer ",
      "rests on the fixed-centroid assumption. Report as a sensitivity analysis ",
      "in Limitations."), max_mdd, SESOI_CINDEX, n_equiv, n_tested)
  } else {
    sprintf(paste0(
      "EQUIVALENCE DOES NOT SURVIVE for all combinations. Largest MDD under the ",
      "full nested bootstrap is %.4f against a SESOI of %.2f; %d/%d view-by-k ",
      "combinations remain equivalent after BH and %d are superior. Once ",
      "discovery-cohort uncertainty is propagated the design cannot BOUND the ",
      "effect for the remaining combinations. The point estimates stay near zero ",
      "-- the change is in precision, not in direction. Abstract, Results and ",
      "Conclusion must be revised to scope the equivalence claim to the ",
      "fixed-classifier estimand."),
      max_mdd, SESOI_CINDEX, n_equiv, n_tested, n_sup)
  }
}
cat("\n=== VERDICT ===\n"); cat(verdict, "\n")

# Relative standard error of an SD estimated from m replicates is
# approximately 1 / sqrt(2(m - 1)). The MDD is proportional to that SD,
# so it inherits the same relative uncertainty.
m_reps <- length(nested_labels)
mdd_rel_se_pct <- if (m_reps > 1) 100 / sqrt(2 * (m_reps - 1)) else NA_real_

writeLines(c(
  "NESTED BOOTSTRAP SENSITIVITY: does the primary equivalence claim",
  "survive propagation of discovery-cohort uncertainty?",
  "======================================================================",
  "",
  "WHY THIS WAS RUN.",
  "Script 14's bootstrap resamples METABRIC patients while holding the",
  "TCGA-derived centroids fixed. That is the correct estimand for",
  "evaluating a pre-trained classifier on a new cohort, and it is the",
  "frame a clinical reader assumes. It does not propagate uncertainty in",
  "the discovery clustering itself. Because this thesis argues that no",
  "null may be claimed without a precision statement, leaving that",
  "assumption's effect on precision unquantified would apply the",
  "standard unevenly. This script quantifies it.",
  "",
  "DESIGN.",
  sprintf("%d replicates requested; %d valid. Each resamples the TCGA cohort",
          N_BOOT_NESTED, length(nested_labels)),
  "with replacement, re-derives the expression-only and fused clusterings",
  "on the resample, rebuilds centroids, re-projects METABRIC, and scores",
  "the paired difference. Two schemes are reported from the same label",
  "sets: DISCOVERY-ONLY (METABRIC fixed) isolates the omitted variance;",
  "FULL NESTED (both resampled) is total uncertainty and is the number",
  "to report.",
  "",
  "SCOPE. Restricted to the two views for which script 07 computed",
  "confidence intervals, fixed before any point estimate was inspected:",
  paste0("  ", paste(VIEWS_NESTED, collapse = ", ")),
  "The question is whether the precision claim holds, not whether a",
  "thirtieth view behaves unusually.",
  "",
  "SESOI. Unchanged at 0.02 C-index. It was fixed before any result was",
  "seen in script 14 and is not revised here. Revising it in the light of",
  "a sensitivity analysis would convert the equivalence framework into",
  "the forking-path exercise it exists to prevent.",
  "",
  "VERDICT.",
  verdict,
  "",
  "LIMITATIONS OF THIS ANALYSIS ITSELF.",
  "1. Feature selection is not resampled. The top-2000 and top-5000",
  "   variance gene sets are held at their full-cohort values.",
  "   Re-selecting features within each replicate would propagate that",
  "   variance too and would widen these intervals further. The MDD",
  "   reported here is a LOWER BOUND on total analytic uncertainty.",
  "2. The MDD is itself an estimate. It is a linear function of the",
  sprintf("   bootstrap standard deviation over %d replicates, whose own",
          length(nested_labels)),
  sprintf("   relative standard error is about %.1f%%. Treat the MDD as",
          mdd_rel_se_pct),
  "   accurate to roughly one significant figure, and do not read a",
  "   narrow margin either side of the SESOI as decisive.",
  "3. n-out-of-n bootstrap with replacement was used, matching scripts",
  "   07 and 14. m-out-of-n subsampling avoids duplicate-sample",
  "   degeneracy in kNN graphs but estimates variance on a different",
  "   scale and would not be comparable without rescaling.",
  "",
  "Generated by scripts/19_nested_bootstrap_sensitivity.R"
), here::here("results/tables/NOTE_nested_bootstrap_sensitivity.txt"))

cat("\nOutputs written:\n",
    " results/tables/table_nested_bootstrap_equivalence.csv\n",
    " results/tables/table_nested_vs_fixed_precision.csv\n",
    " results/tables/NOTE_nested_bootstrap_sensitivity.txt\n",
    " results/objects/bootstrap_nested_discovery.rds\n",
    " results/objects/nested_rep_*.rds  (per-replicate cache)\n")

log_session_info(script_name, key_packages = c("survival", "SNFtool", "org.Hs.eg.db"))
close_logger(log_con, script_name)
