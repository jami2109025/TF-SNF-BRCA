# ============================================================
# 19b_nested_bootstrap_resume_stageC.R
# RESUME script for 19_nested_bootstrap_sensitivity.R
# ============================================================
# WHY THIS SCRIPT EXISTS
# ------------------------------------------------------------
# The 2026-08-17 run of script 19 completed Stage A (100/100 replicates,
# ~6.5h) and Stage B (METABRIC projection + scoring, "Scored replicates:
# 100") successfully -- both are fully cached on disk. It then crashed
# in STAGE C on a bare select() call that resolved to
# AnnotationDbi::select() (an S4 generic pulled in by library(org.Hs.eg.db)
# in stage B) instead of dplyr::select(), because AnnotationDbi has no
# S4 method for a tbl_df:
#
#   Error: unable to find an inherited method for function 'select'
#   for signature 'x = "tbl_df"'
#
# That line has been fixed in 19_nested_bootstrap_sensitivity.R
# (dplyr::select(...) now). But re-sourcing the FULL script would also
# re-run Stage A and Stage B -- both already done and cached -- costing
# another 6.5-9 hours for nothing.
#
# This script starts from the two cache artefacts Stage A/B already
# wrote:
#   - results/objects/nested_rep_%04d.rds   (100 per-replicate label caches)
#   - results/objects/bootstrap_nested_discovery.rds  (boot_nested, the
#     scored paired differences -- written at line 802 of script 19,
#     confirmed on disk since "Scored replicates: 100" printed in the
#     crashed log AFTER the saveRDS() call that produces it)
#
# and reconstructs exactly the two objects Stage C needs (nested_labels,
# boot_nested) from them -- no SNF, no distance matrices, no cluster,
# no METABRIC re-projection. Everything from "STAGE C" in the original
# script onward is then run unchanged (with the select() fix applied).
#
# DELIBERATELY NOT LOADED: org.Hs.eg.db. Nothing downstream of Stage B
# needs it, and loading it is exactly what caused the crash being
# repaired here -- so it is left out rather than fixed-and-reloaded.
#
# HOW TO RUN
# ------------------------------------------------------------
# Place this file in the same scripts/ folder as
# 19_nested_bootstrap_sensitivity.R (so here::here() resolves to the
# same project root, D:/TF-SNF-BRCA) and source it the same way:
#
#   source("D:/TF-SNF-BRCA/scripts/19b_nested_bootstrap_resume_stageC.R")
#
# It writes to the SAME output paths as script 19, so anything
# downstream (report generation, etc.) picks the results up exactly as
# if script 19 had completed on its own. Expected wall-clock: seconds
# to low minutes.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})
source(here::here("R", "utils.R"))
source(here::here("R", "utils_benchmark.R"))

script_name <- "19b_nested_bootstrap_resume_stageC"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()   # matches script 19's discipline; nothing below
                       # consumes fresh randomness -- boot_nested's
                       # values are already fixed on disk.

fail_clean <- function(...) {
  msg <- paste0(...)
  try(close_logger(log_con, script_name), silent = TRUE)
  stop(msg, call. = FALSE)
}

# ------------------------------------------------------------------
# CONFIG -- must match 19_nested_bootstrap_sensitivity.R exactly.
# These are not re-derived from anything; they are the same constants
# the crashed run used, needed here only to label output and to locate
# the same cache files.
# ------------------------------------------------------------------
N_BOOT_NESTED <- 100
K_VALUES      <- c(2, 3)
SESOI_CINDEX  <- 0.02
VIEWS_NESTED  <- c("dorothea_AC__viper", "collectri__viper")

# ------------------------------------------------------------------
# RECONSTRUCT nested_labels FROM THE PER-REPLICATE CACHE
# (cheap: 100 small RDS files holding cluster labels + resample
# indices, not the distance matrices or SNF output -- this is the same
# rep_path()/read loop as script 19 lines 434-436 and 574-589, run
# stand-alone here.)
# ------------------------------------------------------------------
rep_path <- function(b) {
  here::here("results/objects", sprintf("nested_rep_%04d.rds", b))
}

all_b <- seq_len(N_BOOT_NESTED)
missing_reps <- all_b[!file.exists(vapply(all_b, rep_path, character(1)))]
if (length(missing_reps)) {
  fail_clean(
    "Expected all ", N_BOOT_NESTED, " per-replicate caches on disk but ",
    length(missing_reps), " are missing (b = ",
    paste(head(missing_reps, 10), collapse = ", "),
    if (length(missing_reps) > 10) ", ..." else "",
    "). Stage A did not actually finish for these -- this resume script ",
    "cannot proceed without re-running them. Fall back to sourcing the ",
    "full 19_nested_bootstrap_sensitivity.R, which will skip everything ",
    "already cached and compute only these.")
}

nested_labels <- lapply(all_b, function(b) readRDS(rep_path(b)))
ok_flag <- vapply(nested_labels, function(x) isTRUE(x$ok), logical(1))
if (any(!ok_flag)) {
  cat("Failed replicates and reasons:\n")
  for (r in nested_labels[!ok_flag]) cat("  b =", r$b, ":", r$reason, "\n")
}
cat(sprintf("Loaded from cache: %d/%d replicates valid.\n",
            sum(ok_flag), length(ok_flag)))
if (sum(ok_flag) < 50) {
  fail_clean("Fewer than 50 valid replicates; tost_equivalence_bootstrap() ",
             "will refuse to run and no equivalence claim would be meaningful.")
}
nested_labels <- nested_labels[ok_flag]

# ------------------------------------------------------------------
# RECONSTRUCT boot_nested FROM THE STAGE-B CACHE
# ------------------------------------------------------------------
boot_nested_path <- here::here("results/objects/bootstrap_nested_discovery.rds")
if (!file.exists(boot_nested_path)) {
  fail_clean("bootstrap_nested_discovery.rds not found at ", boot_nested_path,
             " -- Stage B did not finish writing it. This resume script ",
             "cannot substitute for Stage B; re-run 19_nested_bootstrap_",
             "sensitivity.R, which will skip Stage A (cached) and redo ",
             "only Stage B (~10-20 min).")
}
boot_nested <- readRDS(boot_nested_path)
cat("Loaded cached boot_nested:", length(unique(boot_nested$b)),
    "scored replicates,", nrow(boot_nested), "rows.\n")

# ==================================================================
# EVERYTHING BELOW IS UNCHANGED FROM 19_nested_bootstrap_sensitivity.R,
# STAGE C ONWARD (its lines ~806-984), with the select() fix applied.
# ==================================================================

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
print(equiv_nested %>% filter(endpoint == "cindex_adjusted_full_nested") %>%
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
  "Generated by scripts/19b_nested_bootstrap_resume_stageC.R",
  "(Stage A/B outputs reused unchanged from the 2026-08-17 run of",
  " scripts/19_nested_bootstrap_sensitivity.R)"
), here::here("results/tables/NOTE_nested_bootstrap_sensitivity.txt"))

cat("\nOutputs written:\n",
    " results/tables/table_nested_bootstrap_equivalence.csv\n",
    " results/tables/table_nested_vs_fixed_precision.csv\n",
    " results/tables/NOTE_nested_bootstrap_sensitivity.txt\n",
    " (results/objects/bootstrap_nested_discovery.rds and nested_rep_*.rds\n",
    "  were reused as-is, not rewritten)\n")

log_session_info(script_name, key_packages = c("here", "tidyverse"))
close_logger(log_con, script_name)
