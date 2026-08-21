# ============================================================
# 00b_verify_blas_swap.R
# Verify an optimised-BLAS installation BEFORE trusting any result
# computed with it. Run once before the swap, once after.
# ============================================================
# WHY VERIFICATION IS NOT OPTIONAL
# ------------------------------------------------------------
# Replacing Rblas.dll changes the numerical library underneath every
# matrix operation in the pipeline. Three things can go wrong, and only
# the first announces itself:
#
#   1. R fails to start        -- obvious, and fixed by rolling back.
#   2. R starts but computes WRONG ANSWERS. This happens if the wrong
#      OpenBLAS variant is installed -- specifically the ILP64 build
#      (64-bit integer indices), which R cannot use. It does not crash;
#      it silently returns garbage for large matrices. Every downstream
#      number would be wrong with no warning.
#   3. R starts, computes correctly, but runs SLOWER than before,
#      because a multithreaded BLAS inside 7 parallel workers
#      oversubscribes the CPU.
#
# This script tests for all three. Do not skip it because R started
# successfully -- case 2 is the dangerous one, and starting successfully
# is exactly what it does.
#
# It also answers the question that matters scientifically: does
# changing BLAS change THIS PIPELINE'S RESULTS? A different BLAS sums
# floating-point numbers in a different order, so tiny differences
# (~1e-12) are EXPECTED and harmless. What must not change is the
# clustering. The final check rebuilds the primary fused network and
# confirms the cluster assignments are identical (ARI = 1).
# ============================================================

cat("=========================================================\n")
cat("BLAS verification\n")
cat("=========================================================\n\n")

cat("R version : ", R.version.string, "\n", sep = "")
cat("R bin dir : ", R.home("bin"), "\n\n", sep = "")

# ------------------------------------------------------------------
# TEST 1. Speed
# ------------------------------------------------------------------
cat("--- TEST 1: matrix multiply speed ---\n")
set.seed(1)
n <- 1095
A <- matrix(rnorm(n * n), n)
el <- system.time(A %*% A)[["elapsed"]]
gflops <- 2 * n^3 / el / 1e9
cat(sprintf("  %d x %d matmul: %.3f s  (%.1f GFLOPS)\n", n, n, el, gflops))
if (gflops < 2) {
  cat("  -> REFERENCE BLAS still in use. The swap has not taken effect.\n")
  cat("     Check that Rblas.dll was actually replaced and R was restarted.\n")
} else {
  cat("  -> OPTIMISED BLAS active.\n")
  cat(sprintf("     Expected SNF time: ~%.0f s (was ~696 s)\n", 696 * (9.26 / el)))
}

# ------------------------------------------------------------------
# TEST 2. Numerical correctness -- the ILP64 trap
#
# An ILP64 OpenBLAS passes small tests and fails on large ones, because
# the index overflow only bites past a certain matrix size. Both sizes
# are therefore tested. If the large case fails while the small one
# passes, the wrong OpenBLAS variant is installed: download the
# standard x64 build, NOT the one whose filename ends in _64.
# ------------------------------------------------------------------
cat("\n--- TEST 2: numerical correctness ---\n")
check_identity <- function(m) {
  set.seed(42)
  X <- matrix(rnorm(m * m), m)
  err <- max(abs(X %*% solve(X) - diag(m)))
  ok <- is.finite(err) && err < 1e-6
  cat(sprintf("  %4d x %4d : max|X X^-1 - I| = %.3e  %s\n",
              m, m, err, if (ok) "PASS" else "*** FAIL ***"))
  ok
}
ok_small <- check_identity(200)
ok_large <- check_identity(1500)

if (!ok_large) {
  cat("\n  *** STOP. The BLAS is returning incorrect results on large\n")
  cat("      matrices. This is the signature of an ILP64 (64-bit integer)\n")
  cat("      OpenBLAS build, which R cannot use. Roll back to\n")
  cat("      Rblas.dll.BACKUP and download the standard x64 build --\n")
  cat("      the one WITHOUT '_64' in the filename.\n")
  cat("      DO NOT run any analysis until this passes.\n")
  stop("BLAS numerical check failed.")
}

# Eigen decomposition exercises LAPACK, which calls through to BLAS --
# spectralClustering() depends on it, so it is checked separately.
set.seed(7)
S <- crossprod(matrix(rnorm(600 * 600), 600))
ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
cat(sprintf("  eigen() on 600x600 SPD matrix: all eigenvalues positive and finite: %s\n",
            all(is.finite(ev)) && all(ev > 0)))

# ------------------------------------------------------------------
# TEST 3. Threading -- does BLAS oversubscribe under parallel workers?
# ------------------------------------------------------------------
cat("\n--- TEST 3: BLAS threading ---\n")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  cat("  RhpcBLASctl installed. BLAS threads reported:",
      tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NA), "\n")
  cat("  Physical cores:", tryCatch(RhpcBLASctl::get_num_cores(), error = function(e) NA), "\n")
  cat("  -> scripts 07 and 14 pin this to 1 inside each worker; that is correct.\n")
} else {
  cat("  RhpcBLASctl NOT installed. Install it before running scripts 07/14:\n")
  cat("      install.packages('RhpcBLASctl')\n")
  cat("  Without it, a multithreaded BLAS inside 7 parallel workers can run\n")
  cat("  SLOWER than the reference BLAS did.\n")
}

# ------------------------------------------------------------------
# TEST 4. Does the pipeline's own result change?
#
# The decisive scientific question. A different BLAS sums in a different
# order, so the fused network will differ in the last few decimal places
# -- that is expected and harmless. What must NOT change is the
# partition. ARI = 1 against the stored labels means the swap is
# scientifically transparent and results computed before and after can
# be reported together.
# ------------------------------------------------------------------
cat("\n--- TEST 4: does the pipeline's clustering change? ---\n")
have <- file.exists(here::here("results/objects/W_fused.rds")) &&
        file.exists(here::here("data/processed/vst_top2000_genes.rds")) &&
        file.exists(here::here("data/processed/tf_activity_viper_AC_primary.rds"))

if (!have) {
  cat("  Skipped: run scripts 05 and 06 first.\n")
} else {
  suppressPackageStartupMessages({ library(here); library(SNFtool); library(mclust) })
  source(here::here("R", "utils.R"))
  source(here::here("R", "utils_benchmark.R"))

  expr <- readRDS(here::here("data/processed/vst_top2000_genes.rds"))
  tfa  <- readRDS(here::here("data/processed/tf_activity_viper_AC_primary.rds"))
  Wref <- readRDS(here::here("results/objects/W_fused.rds"))

  common <- Reduce(intersect, list(colnames(expr), colnames(tfa), colnames(Wref)))
  expr <- expr[, common]; tfa <- tfa[, common]; Wref <- Wref[common, common]

  sc <- function(m) { s <- scale(t(m)); s[is.na(s)] <- 0; s }
  W_e <- affinity_from_dist(dist2_cached(sc(expr)), K = 20, sigma = 0.5)
  W_t <- affinity_from_dist(dist2_cached(sc(tfa)),  K = 20, sigma = 0.5)
  Wnew <- SNFtool::SNF(list(W_e, W_t), K = 20, t = 20)

  dev <- max(abs(Wnew - Wref))
  cat(sprintf("  max |rebuilt - stored| in fused network: %.3e\n", dev))
  cat("    (a value around 1e-12 is EXPECTED after a BLAS change -- different\n")
  cat("     summation order, same mathematics. Exact 0 would mean the BLAS\n")
  cat("     did not actually change.)\n")

  for (k in c(2, 3)) {
    a <- SNFtool::spectralClustering(Wref, K = k)
    b <- SNFtool::spectralClustering(Wnew, K = k)
    ari <- mclust::adjustedRandIndex(a, b)
    cat(sprintf("  k = %d : ARI(stored vs rebuilt) = %.6f  %s\n", k, ari,
                if (ari > 0.999) "IDENTICAL PARTITION" else "*** PARTITION CHANGED ***"))
    if (ari <= 0.999) {
      cat("      -> The clustering is not stable to the BLAS change. Do not mix\n")
      cat("         results computed before and after the swap. Re-run script 06\n")
      cat("         under the new BLAS so the whole pipeline is internally\n")
      cat("         consistent, and say which BLAS was used in the methods.\n")
    }
  }
}

cat("\n=========================================================\n")
cat("If TEST 2 passed and TEST 4 shows ARI = 1, the swap is safe.\n")
cat("Record the BLAS in the methods section either way -- it is part\n")
cat("of the computational environment and belongs in the archive\n")
cat("produced by 99_archive_external_resources.R.\n")
cat("=========================================================\n")
