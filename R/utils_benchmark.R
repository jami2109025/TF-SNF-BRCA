# ============================================================
# R/utils_benchmark.R
# Helpers for the ADDED-VALUE benchmarking framework.
#
# This file is DELIBERATELY SEPARATE from R/utils.R so that the
# previously-audited helper set stays byte-identical and its audit
# trail remains valid. Source BOTH files in any script that needs
# these functions:
#
#     source(here::here("R", "utils.R"))
#     source(here::here("R", "utils_benchmark.R"))
#
# ------------------------------------------------------------
# WHY THIS FILE EXISTS
# ------------------------------------------------------------
# The original pipeline asked "does adding a TF-activity view improve
# clustering?", answered it with ONE regulon x ONE inference method x
# ONE agreement metric (ARI vs PAM50), obtained a confidence interval
# spanning zero, and had no way to distinguish the two very different
# conclusions that a null result permits:
#
#   (a) "TF activity genuinely adds nothing"                 <- a finding
#   (b) "this particular configuration/metric had no power"  <- a non-finding
#
# Only (a) is publishable. Separating (a) from (b) requires three
# things that the original design did not have, all provided here:
#
#   1. EQUIVALENCE TESTING (tost_equivalence_bootstrap). A confidence
#      interval containing zero is NOT evidence of no effect. It only
#      becomes evidence of no effect once the interval is also shown to
#      exclude every effect large enough to matter. That requires a
#      PRE-SPECIFIED smallest effect size of interest (SESOI) and a
#      two-one-sided-test procedure against it.
#
#   2. PAIRED RESAMPLING (make_resample_index). Comparing view A and
#      view B on independently drawn bootstrap samples inflates the
#      variance of their difference enormously. Every view must be
#      evaluated on the SAME resample so the difference is a paired
#      contrast. Achieved here by making the resample index a pure
#      deterministic function of the replicate number b.
#
#   3. A CACHED DISTANCE PATH (dist2_cached / affinity_from_dist).
#      Pairwise sample distance does not depend on which samples were
#      drawn -- d(i,j) is a property of the pair. So the O(n^2 * p)
#      distance computation can be done ONCE on the full cohort and
#      every resample can index into it, D[idx, idx]. Only the
#      O(n^2 log n) affinity/kNN step must be recomputed. This is what
#      makes a multi-view bootstrap feasible at all: the original
#      single-view bootstrap took 22.8 hours.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

# ------------------------------------------------------------------
# 0. PRE-SPECIFIED ANALYSIS CONSTANTS
#
# These are declared here, in one place, rather than inline at the
# point of use, precisely so that they are visibly PRE-specified
# rather than chosen after seeing the results. Any change to these
# values is a change to the analysis plan and should be recorded as
# such (git history is the audit trail).
# ------------------------------------------------------------------

#' Smallest effect size of interest for the added-value contrast,
#' on the ARI scale.
#'
#' JUSTIFICATION (state this in the manuscript, do not leave it
#' unargued -- an unjustified SESOI is the single easiest thing for a
#' reviewer to attack in an equivalence-testing paper):
#'
#'   - The ARI difference between this pipeline's own k=2 and k=3
#'     solutions against PAM50 is 0.411 - 0.326 = 0.085. An integration
#'     benefit smaller than the effect of changing k by one is not a
#'     meaningful reason to prefer the more complex model.
#'   - Multi-omic integration papers that report a benefit typically
#'     report ARI/NMI gains in the 0.05-0.20 range; 0.05 sits at the
#'     bottom of that range, making this a CONSERVATIVE (i.e. hard to
#'     pass) equivalence bound.
#'   - Below ~0.05 ARI the added view cannot change the clinical or
#'     biological interpretation of any cluster, because it corresponds
#'     to a relabelling of a few tens of samples out of ~1100.
#'
#' A conservative bound matters here: passing equivalence against a
#' SMALL bound is a STRONG claim. Passing against a large bound is
#' nearly vacuous.
SESOI_ARI <- 0.05

#' Equivalence-test alpha. The 90% CI rule below is the standard
#' CI-based form of TOST at alpha = 0.05 (two one-sided tests at 0.05
#' each correspond to a 1 - 2*alpha = 90% interval).
EQUIV_ALPHA <- 0.05

#' Numerical guard for degenerate affinity construction under
#' with-replacement resampling (see affinity_from_dist).
AFFINITY_EPS <- 1e-8

# ------------------------------------------------------------------
# 1. Cached distance / affinity path
# ------------------------------------------------------------------

#' Squared Euclidean distance matrix, computed ONCE for a full cohort.
#'
#' Deliberately matches SNFtool::dist2()'s output convention (SQUARED
#' Euclidean distance) so that cached distances are drop-in compatible
#' with SNFtool::affinityMatrix(). Verified against SNFtool::dist2() by
#' verify_distance_cache() below -- do not skip that verification, the
#' entire speedup is invalid if the conventions diverge.
#'
#' @param mat samples-in-ROWS, features-in-columns (i.e. already
#'   transposed and scaled, exactly what the pipeline feeds to dist2).
dist2_cached <- function(mat) {
  mat <- as.matrix(mat)
  SNFtool::dist2(mat, mat)
}

#' Hard verification that the cached distance matrix is identical to
#' what the original (uncached) code path would have produced.
#'
#' This is the correctness guarantee for the entire caching
#' optimisation: if D[idx, idx] equals dist2(mat[idx, ], mat[idx, ])
#' then every resampled network built from the cache is bit-identical
#' to one built the slow way. Run it on a small random subset -- it is
#' cheap, and it converts "this optimisation should be safe" into
#' "this optimisation was checked at runtime".
verify_distance_cache <- function(mat, D_full, n_check = 60, tol = 1e-8, label = "") {
  mat <- as.matrix(mat)
  set_pipeline_seed(offset = 9001)
  idx <- sample(seq_len(nrow(mat)), min(n_check, nrow(mat)))
  D_direct <- SNFtool::dist2(mat[idx, , drop = FALSE], mat[idx, , drop = FALSE])
  D_cached <- D_full[idx, idx, drop = FALSE]
  max_abs_diff <- max(abs(D_direct - D_cached))
  ok <- max_abs_diff < tol
  cat("[cache check", label, "] max |direct - cached| =", format(max_abs_diff, digits = 3),
      if (ok) " -- PASS\n" else " -- FAIL\n")
  if (!ok) {
    stop("verify_distance_cache(): cached distance matrix does NOT reproduce ",
         "SNFtool::dist2() output for label '", label, "'. The caching ",
         "optimisation is invalid under the current SNFtool version -- do not ",
         "trust any bootstrap result from this run.")
  }
  invisible(TRUE)
}

#' Affinity matrix from a PRE-COMPUTED distance matrix.
#'
#' DEGENERACY GUARD: under with-replacement bootstrap resampling a
#' sample can be drawn more than once, producing exact zero distances
#' between duplicated copies. SNFtool::affinityMatrix() normalises each
#' row by the mean distance to its K nearest neighbours; if enough of
#' those neighbours are zero-distance duplicates, that mean is 0 and
#' the affinity is NaN/Inf. The original bootstrap had this exposure
#' and no guard -- a replicate silently producing NaN would propagate
#' into spectralClustering() and could be dropped or, worse, produce
#' a meaningless partition that still returned a finite ARI.
#'
#' Here we (i) add a tiny epsilon to exact-zero off-diagonal
#' distances, breaking only the degeneracy and nothing else, and
#' (ii) hard-check the result is finite.
affinity_from_dist <- function(D, K = 20, sigma = 0.5) {
  D <- as.matrix(D)
  off_diag_zero <- (D == 0)
  diag(off_diag_zero) <- FALSE
  if (any(off_diag_zero)) {
    D[off_diag_zero] <- AFFINITY_EPS
  }
  W <- SNFtool::affinityMatrix(D, K = K, sigma = sigma)
  if (!all(is.finite(W))) {
    stop("affinity_from_dist(): non-finite affinity values produced even after ",
         "the zero-distance degeneracy guard. Inspect the input distance matrix.")
  }
  W
}

#' Deterministic resample index for replicate b.
#'
#' CRITICAL for the paired design: because the index depends only on b
#' (and the pipeline seed), EVERY view evaluated at replicate b sees
#' the SAME resampled patients. The added-value difference is therefore
#' a paired contrast within replicate, and the bootstrap distribution
#' of that difference reflects only the view effect, not resampling
#' noise shared by both arms.
#'
#' @param scheme "bootstrap" = n draws with replacement (classic,
#'   matches the original pipeline's design so results stay
#'   comparable). "subsample" = m-out-of-n WITHOUT replacement, which
#'   cannot create duplicate-sample degeneracy at all and is the
#'   cleaner choice for kNN-graph methods; reported as a sensitivity
#'   analysis.
make_resample_index <- function(b, n, scheme = c("bootstrap", "subsample"),
                                subsample_frac = 0.8, seed_base = 2000L) {
  scheme <- match.arg(scheme)
  set.seed(PIPELINE_SEED + seed_base + b)
  if (scheme == "bootstrap") {
    sample(seq_len(n), n, replace = TRUE)
  } else {
    sample(seq_len(n), round(subsample_frac * n), replace = FALSE)
  }
}

# ------------------------------------------------------------------
# 2. Regulon handling
# ------------------------------------------------------------------

#' Normalise any regulon/network source to a common schema.
#'
#' DoRothEA uses columns (tf, confidence, target, mor, likelihood);
#' CollecTRI uses (source, target, mor) with NO likelihood column.
#' Downstream code should never have to branch on which source it got,
#' so everything is coerced to (source, target, mor, likelihood) here,
#' with likelihood defaulting to 1 (uniform weighting) when the source
#' does not supply interaction confidence.
normalise_network <- function(net) {
  net <- as_tibble(net)
  if ("tf" %in% colnames(net) && !"source" %in% colnames(net)) {
    net <- dplyr::rename(net, source = tf)
  }
  if (!all(c("source", "target", "mor") %in% colnames(net))) {
    stop("normalise_network(): expected columns source/tf, target, mor. Got: ",
         paste(colnames(net), collapse = ", "))
  }
  if (!"likelihood" %in% colnames(net)) {
    net$likelihood <- 1
  }
  net %>%
    dplyr::select(source, target, mor, likelihood) %>%
    dplyr::filter(!is.na(source), !is.na(target)) %>%
    dplyr::distinct(source, target, .keep_all = TRUE)
}

#' Degree-preserving edge permutation: the NETWORK-SIDE negative control.
#'
#' This is a DIFFERENT and strictly stronger control than the sample
#' permutation already in the pipeline, and both are needed because
#' they falsify different things:
#'
#'   - SAMPLE permutation (original pipeline): shuffles which patient
#'     each TF-activity profile belongs to. Destroys the TF-to-patient
#'     link. Answers "is fusion responding to real per-patient TF
#'     structure, or would any noise view do?"
#'
#'   - EDGE permutation (here): keeps the patients, keeps each TF's
#'     out-degree and the mor/likelihood distribution, but reassigns
#'     WHICH GENES each TF regulates. Answers the much sharper
#'     question "is the benefit (if any) coming from CURATED REGULATORY
#'     BIOLOGY, or merely from the dimensionality reduction of
#'     averaging a few hundred random gene sets?"
#'
#' The second question is the one a sceptical reviewer will actually
#' ask, because a randomised-regulon control that performs as well as
#' the real regulon would show that "TF activity" is functioning as a
#' generic smoothing operation rather than as regulatory inference.
permute_network_edges <- function(net, seed_offset = 4200L) {
  net <- normalise_network(net)
  set_pipeline_seed(offset = seed_offset)
  all_targets <- net$target
  net$target <- sample(all_targets)          # global target relabelling
  net %>% dplyr::distinct(source, target, .keep_all = TRUE)
}

#' Filter a network to targets measured in the matrix and to TFs with
#' at least `min_targets` measured targets. Returns network + size table.
filter_network_to_matrix <- function(net, mat, min_targets = 10) {
  net <- normalise_network(net) %>% dplyr::filter(target %in% rownames(mat))
  sizes <- net %>% dplyr::count(source, name = "n_targets") %>%
    dplyr::filter(n_targets >= min_targets)
  list(
    network = net %>% dplyr::filter(source %in% sizes$source),
    sizes   = sizes
  )
}

# ------------------------------------------------------------------
# 3. Evaluation endpoints
#
# The original pipeline used exactly ONE endpoint: ARI against PAM50.
# That single choice is the deepest design flaw in the original
# analysis, for a reason that is easy to miss:
#
#   PAM50 labels are themselves computed FROM gene expression. An
#   expression-only clustering is therefore near-optimal for that
#   endpoint BY CONSTRUCTION. The benchmark was built so that the
#   baseline arm is playing at home and the TF arm can essentially
#   only lose. A null result on a rigged endpoint is uninformative.
#
# The fix is to evaluate on endpoints that are NOT downstream of the
# same expression matrix. Below, in increasing order of independence:
#
#   ari_vs_reference / nmi_vs_reference -- expression-derived label
#       agreement. RETAINED for continuity with the original result and
#       with the wider literature, but explicitly flagged as circular.
#
#   cox_lrt_endpoint -- prognostic separation. Survival outcome is
#       NOT derived from the expression matrix, so this endpoint does
#       not structurally favour either arm. This is the endpoint that
#       can actually be won or lost on merit.
#
#   concordance_endpoint -- discrimination (Harrell's C) on the same
#       outcome, on an interpretable 0.5-1 scale, which is what a
#       clinical reader will want to see.
# ------------------------------------------------------------------

#' Align two label vectors, dropping any element missing in either.
#'
#' WHY THIS EXISTS. Under with-replacement bootstrap resampling a patient
#' can be drawn several times, so the affinity graph contains
#' near-duplicate nodes and spectralClustering() occasionally returns NA
#' for a handful of them. That is rare and harmless to a summary
#' statistic -- but aricode::NMI() responds to a single NA by throwing
#' "NA are not supported", which inside parLapply kills every worker and
#' takes the entire run with it. A secondary endpoint must never be able
#' to destroy a multi-day job, so every endpoint below now filters
#' missing values explicitly and returns NA rather than signalling.
.clean_pair <- function(cl, reference, keep) {
  if (is.null(keep)) keep <- rep(TRUE, length(cl))
  keep <- keep & !is.na(keep)
  a <- cl[keep]; b <- reference[keep]
  ok <- !is.na(a) & !is.na(b)
  list(a = a[ok], b = b[ok], n = sum(ok))
}

ari_vs_reference <- function(cl, reference, keep) {
  p <- .clean_pair(cl, reference, keep)
  if (p$n <= 10) return(NA_real_)
  tryCatch(mclust::adjustedRandIndex(p$a, p$b), error = function(e) NA_real_)
}

nmi_vs_reference <- function(cl, reference, keep) {
  if (!requireNamespace("aricode", quietly = TRUE)) return(NA_real_)
  p <- .clean_pair(cl, reference, keep)
  if (p$n <= 10) return(NA_real_)
  # Integer factor codes, not character: aricode coerces character input
  # internally and emits "NAs introduced by coercion" on non-numeric
  # labels such as "LumA". Passing codes avoids that path entirely.
  a <- as.integer(factor(as.character(p$a)))
  b <- as.integer(factor(as.character(p$b)))
  if (length(unique(a)) < 2 || length(unique(b)) < 2) return(NA_real_)
  tryCatch(aricode::NMI(a, b), error = function(e) NA_real_)
}

#' Likelihood-ratio chi-square for the cluster term in a Cox model.
#'
#' Scale note: this is a chi-square on (n_clusters - 1) df, so it is
#' comparable ACROSS METHODS at a fixed k (all arms produce k groups)
#' but NOT across different k. Every comparison in this pipeline is
#' made within a fixed k, so that is sufficient.
cox_lrt_endpoint <- function(cl, time, status) {
  ok <- is.finite(time) & time > 0 & !is.na(status) & !is.na(cl)
  ok[is.na(ok)] <- FALSE
  if (sum(ok) < 30 || length(unique(cl[ok])) < 2) return(NA_real_)
  fit <- try(
    survival::coxph(survival::Surv(time[ok], status[ok]) ~ factor(cl[ok])),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(NA_real_)
  unname(summary(fit)$logtest["test"])
}

#' Harrell's C for a cluster-only Cox model (optionally adjusted).
concordance_endpoint <- function(cl, time, status, covariates = NULL) {
  ok <- is.finite(time) & time > 0 & !is.na(status) & !is.na(cl)
  ok[is.na(ok)] <- FALSE
  if (!is.null(covariates)) ok <- ok & stats::complete.cases(covariates)
  if (sum(ok) < 30 || length(unique(cl[ok])) < 2) return(NA_real_)
  df <- data.frame(time = time[ok], status = status[ok], cl = factor(cl[ok]))
  form <- "survival::Surv(time, status) ~ cl"
  if (!is.null(covariates)) {
    cov_ok <- covariates[ok, , drop = FALSE]
    df <- cbind(df, cov_ok)
    form <- paste(form, "+", paste(colnames(cov_ok), collapse = " + "))
  }
  fit <- try(survival::coxph(stats::as.formula(form), data = df), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  unname(summary(fit)$concordance[1])
}

# ------------------------------------------------------------------
# 4. Equivalence testing -- the statistical core of a defensible
#    negative result
# ------------------------------------------------------------------

#' Bootstrap two-one-sided-tests (TOST) for equivalence to zero.
#'
#' A null hypothesis significance test can only ever FAIL TO REJECT the
#' null; it can never accept it. "The 95% CI included zero" is
#' therefore compatible with both "no effect" and "underpowered study",
#' and a reviewer is entitled to assume the latter. Equivalence testing
#' inverts the logic: the null becomes "the effect is at least as large
#' as the SESOI", and REJECTING that null is positive evidence that any
#' true effect is too small to matter.
#'
#' Two decision rules are reported because they are complementary:
#'
#'   (1) CI rule: if the (1 - 2*alpha) interval -- 90% at alpha = 0.05 --
#'       lies entirely inside (-SESOI, +SESOI), equivalence is declared.
#'       This is the standard, transparent, referee-friendly form.
#'
#'   (2) Percentile p-values: p_upper = Pr(delta* >= SESOI) and
#'       p_lower = Pr(delta* <= -SESOI) under the bootstrap
#'       distribution; the TOST p-value is their maximum.
#'
#' The FOUR possible outcomes are enumerated explicitly in the returned
#' `conclusion` field, because the two that are usually conflated --
#' "equivalent" and "inconclusive" -- are exactly the distinction this
#' whole exercise exists to make.
tost_equivalence_bootstrap <- function(diffs, sesoi = SESOI_ARI, alpha = EQUIV_ALPHA,
                                       label = "") {
  diffs <- diffs[is.finite(diffs)]
  n_valid <- length(diffs)
  if (n_valid < 50) {
    return(tibble(
      label = label, n_boot_valid = n_valid,
      point_estimate = NA_real_, ci90_lower = NA_real_, ci90_upper = NA_real_,
      ci95_lower = NA_real_, ci95_upper = NA_real_,
      sesoi = sesoi, p_tost = NA_real_, p_two_sided = NA_real_,
      equivalent = NA, superior = NA,
      conclusion = "insufficient valid bootstrap replicates"
    ))
  }

  point <- mean(diffs)
  ci90  <- stats::quantile(diffs, c(alpha, 1 - alpha), na.rm = TRUE)
  ci95  <- stats::quantile(diffs, c(0.025, 0.975), na.rm = TRUE)

  p_upper <- mean(diffs >= sesoi)
  p_lower <- mean(diffs <= -sesoi)
  p_tost  <- max(p_upper, p_lower)

  # Two-sided bootstrap p-value against zero (for the superiority arm).
  p_two_sided <- 2 * min(mean(diffs <= 0), mean(diffs >= 0))
  p_two_sided <- min(1, p_two_sided)

  equivalent <- (ci90[1] > -sesoi) && (ci90[2] < sesoi)
  superior   <- ci95[1] > 0

  conclusion <- dplyr::case_when(
    superior &  equivalent ~ "significant but trivially small (statistically > 0, practically within the SESOI)",
    superior & !equivalent ~ "SUPERIOR: added value is significantly greater than zero",
   !superior &  equivalent ~ "EQUIVALENT TO ZERO: any true added value is smaller than the pre-specified SESOI",
    TRUE                   ~ "INCONCLUSIVE: cannot exclude zero, and cannot exclude a meaningful effect either -- underpowered"
  )

  tibble(
    label = label, n_boot_valid = n_valid,
    point_estimate = point,
    ci90_lower = unname(ci90[1]), ci90_upper = unname(ci90[2]),
    ci95_lower = unname(ci95[1]), ci95_upper = unname(ci95[2]),
    sesoi = sesoi, p_tost = p_tost, p_two_sided = p_two_sided,
    equivalent = equivalent, superior = superior,
    conclusion = conclusion
  )
}

#' Multiplicity control across a grid of views.
#'
#' Testing ~15 regulon x method views for added value and reporting the
#' best one is a garden-of-forking-paths problem: with 15 independent
#' tests at alpha = 0.05 the probability of at least one spurious
#' "significant" view is 1 - 0.95^15 = 54%. Any view-level superiority
#' claim must therefore be BH-adjusted across the whole grid, and the
#' adjusted column is what gets reported in the manuscript.
#'
#' Note the asymmetry, which is a feature and worth stating in the
#' paper: multiplicity inflates FALSE POSITIVES, so it makes a
#' superiority claim harder while making the overall NEGATIVE
#' conclusion more conservative, not less. If no view survives
#' adjustment, that is a stronger negative result, not a weaker one.
adjust_grid_multiplicity <- function(equiv_table) {
  equiv_table %>%
    dplyr::mutate(
      p_two_sided_BH = stats::p.adjust(p_two_sided, method = "BH"),
      p_tost_BH      = stats::p.adjust(p_tost,      method = "BH"),
      superior_after_BH  = p_two_sided_BH < 0.05 & point_estimate > 0,
      equivalent_after_BH = p_tost_BH < 0.05
    )
}

#' Minimum detectable difference at the observed bootstrap precision.
#'
#' Answers, in advance of any reviewer asking, "was this study even
#' capable of detecting the effect it failed to find?". If the MDD is
#' larger than the SESOI, the design was underpowered and no negative
#' claim should be made regardless of what the CI looks like.
minimum_detectable_difference <- function(diffs, power = 0.8, alpha = 0.05) {
  diffs <- diffs[is.finite(diffs)]
  if (length(diffs) < 50) return(NA_real_)
  se <- stats::sd(diffs)
  (stats::qnorm(1 - alpha / 2) + stats::qnorm(power)) * se
}

# ------------------------------------------------------------------
# 5. BLAS thread control -- required whenever an optimised BLAS is
#    combined with parallel workers
# ------------------------------------------------------------------

#' Force single-threaded BLAS inside a parallel worker.
#'
#' WHY THIS IS NECESSARY, AND WHY OMITTING IT CAN MAKE THINGS SLOWER:
#'
#' An optimised BLAS (OpenBLAS, MKL) multithreads matrix multiplication
#' by default, typically across every physical core. If 7 parallel R
#' workers each launch an 8-thread BLAS on an 8-core machine, the
#' machine is asked to run 56 compute threads. They do not share
#' nicely: the operating system context-switches between them, cache
#' lines are evicted constantly, and total throughput can fall BELOW
#' the single-threaded reference-BLAS baseline. This is the classic
#' oversubscription trap, and it is easy to mistake for "the optimised
#' BLAS did not help".
#'
#' The correct arrangement for this pipeline is parallelism at the
#' REPLICATE level (many independent SNF runs) with BLAS pinned to one
#' thread inside each worker. Replicate-level parallelism scales almost
#' perfectly; BLAS-level parallelism on a 1095x1095 matrix does not.
#' Single-threaded OpenBLAS is still roughly 20-30x faster than
#' reference BLAS, so essentially all of the benefit is retained.
#'
#' Safe to call unconditionally: it is a no-op when RhpcBLASctl is not
#' installed or when the BLAS does not support runtime thread control.
set_blas_threads_single <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
  }
  # Environment fallbacks for BLAS builds without a runtime API.
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1")
  invisible(TRUE)
}

#' Report which BLAS is in use, by measuring it.
#'
#' There is no portable way to ask R which BLAS it links against, so
#' this benchmarks a matrix multiply and infers from the throughput.
#' Printed at the top of long jobs, it makes the difference between a
#' 4-hour run and a 40-hour run visible before the run starts rather
#' than after.
report_blas_speed <- function(n = 1000) {
  set_pipeline_seed(offset = 9100)
  A <- matrix(rnorm(n * n), n)
  el <- system.time(A %*% A)[["elapsed"]]
  gflops <- 2 * n^3 / el / 1e9
  cat(sprintf("BLAS check: %d x %d matmul in %.2f s (%.1f GFLOPS) -- %s\n",
              n, n, el, gflops,
              if (gflops < 2) "REFERENCE BLAS (unoptimised); an optimised BLAS would cut SNF time by 20-50x"
              else "optimised BLAS detected"))
  invisible(gflops)
}

cat("R/utils_benchmark.R loaded. SESOI (ARI) =", SESOI_ARI,
    "| equivalence alpha =", EQUIV_ALPHA, "\n")
