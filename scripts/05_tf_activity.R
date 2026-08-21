# ============================================================
# 05_tf_activity.R  (REVISED: multi-regulon x multi-method grid)
# Transcription factor activity inference
# ============================================================
# WHAT CHANGED AND WHY
# ------------------------------------------------------------
# The previous revision inferred TF activity ONE way: DoRothEA A-C
# regulons scored with VIPER. Script 07 then found that this single
# view added nothing to clustering (added-value 95% CI spanned zero),
# and the thesis had no defence against the obvious objection:
#
#     "You did not find a benefit because you configured TF inference
#      badly. Try CollecTRI. Try a linear model instead of VIPER."
#
# That objection is fatal to a single-configuration negative result and
# unanswerable without new computation. This revision answers it by
# construction: TF activity is now inferred across a GRID of
#
#     {4 regulon sources} x {up to 5 inference methods}
#
# and every downstream benchmark (script 07) is run over the whole
# grid. The negative result then reads "no configuration in a
# systematically varied grid added value", which is a defensible
# scientific claim, rather than "our one configuration did not work",
# which is not.
#
# REGULON SOURCES
#   dorothea_AC   DoRothEA confidence A-C. The original primary view;
#                 retained so every previously reported number can
#                 still be reproduced exactly.
#   dorothea_AB   DoRothEA A-B only. Fewer, better-evidenced edges.
#                 Tests whether the A-C result was diluted by
#                 low-confidence interactions.
#   collectri     CollecTRI. The curated successor resource to
#                 DoRothEA from the same group, and the current
#                 recommended default in decoupleR. Using DoRothEA
#                 ALONE in 2026 is the single most likely reviewer
#                 objection to this manuscript; including CollecTRI
#                 removes it.
#   *_edgeperm    Degree-preserving edge-permuted controls. See
#                 permute_network_edges() in R/utils_benchmark.R for
#                 why this control is scientifically load-bearing:
#                 it separates "curated regulatory biology helps" from
#                 "averaging a few hundred arbitrary gene sets helps".
#
# INFERENCE METHODS (all via decoupleR, on the identical signature)
#   viper   Rank-based enrichment; the original method. Retained.
#   ulm     Univariate linear model. decoupleR's benchmarked
#           top performer for TF activity in most settings.
#   mlm     Multivariate linear model; fits all TFs jointly, so it
#           accounts for shared targets between regulons rather than
#           scoring each TF in isolation.
#   norm_wsum  Permutation-normalised weighted sum.
#   consensus  Rank-based consensus across the above.
#
# METHODOLOGICAL POINT CARRIED FORWARD FROM THE PREVIOUS REVISION:
# VIPER (and every other method here) is run on a per-gene,
# across-sample Z-SCORED expression SIGNATURE, not on raw VST values.
# This is required by the rank/enrichment semantics of aREA and is
# equally appropriate for the linear-model methods. The unscaled VST
# matrix remains the input everywhere else in the pipeline.
#
# RUNTIME: dominated by VIPER, whose cost scales with the number of
# regulons. Measured on this cohort: DoRothEA A-C (263 TFs) 26.5 min,
# DoRothEA A-B (111 TFs) 4.4 min, ULM under 1.5 min in both cases.
# CollecTRI carries ~1200 TFs before size filtering, so its VIPER call
# is the single most expensive cell in the grid. With 3 regulon sources
# x {real, edge-permuted} x {viper, ulm}, budget roughly 4 hours.
#
# NOTE ON BLAS: VIPER and SNF are both dominated by dense matrix
# multiplication. R's default reference BLAS on Windows runs at a few
# hundred MFLOPS; an optimised BLAS (OpenBLAS/MKL) is 20-50x faster and
# shortens this script and script 07 proportionally. Check with:
#     A <- matrix(rnorm(1095*1095), 1095); system.time(A %*% A)
# Over ~1.5 s indicates reference BLAS.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(dorothea)
  library(decoupleR)
  library(org.Hs.eg.db)
  library(pheatmap)
})
source(here::here("R", "utils.R"))
source(here::here("R", "utils_benchmark.R"))

script_name <- "05_tf_activity"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

# ------------------------------------------------------------------
# CONFIG -- pre-specified, edit deliberately
# ------------------------------------------------------------------
MIN_REGULON_SIZE   <- 10      # min measured targets for a TF to be scored
WSUM_PERMUTATIONS  <- 100     # norm_wsum permutations; 0 disables norm_wsum
INCLUDE_EDGEPERM   <- TRUE    # degree-preserving edge-permuted controls
METHODS_REQUESTED  <- c("viper", "ulm", "mlm", "norm_wsum")
COLLECTRI_CACHE    <- here::here("data/processed/network_collectri.rds")

# ------------------------------------------------------------------
# Load FULL filtered, VST-normalized expression matrix
# ------------------------------------------------------------------
vst_all         <- readRDS(here::here("data/processed/vst_matrix_all_genes.rds"))
sample_metadata <- readRDS(here::here("data/processed/sample_metadata_matched.rds"))
stopifnot(identical(colnames(vst_all), rownames(sample_metadata)))
cat("Full filtered expression matrix:", nrow(vst_all), "genes x", ncol(vst_all), "samples\n")

# ------------------------------------------------------------------
# Ensembl -> SYMBOL, dedup by highest variance IN THE VST MATRIX, so
# the SAME representative Ensembl ID is chosen here and in script 11.
# ------------------------------------------------------------------
mapped <- map_to_symbol_dedup(vst_all, variance_reference = vst_all, org_db = org.Hs.eg.db)
expr_symbol_mat <- mapped$matrix
symbol_map      <- mapped$map
cat("Expression matrix after SYMBOL mapping/dedup:", nrow(expr_symbol_mat), "genes\n")
write_csv(symbol_map, here::here("results/tables/table_ensembl_symbol_map_used.csv"))

# ------------------------------------------------------------------
# Per-sample gene expression SIGNATURE (per-gene z-score across
# samples). Built ONCE and shared by every regulon x method cell, so
# that grid cells differ ONLY in regulon and scoring method -- never
# in their input data. Any difference observed across the grid is
# therefore attributable to the varied factor and nothing else.
# ------------------------------------------------------------------
expr_signature <- t(scale(t(expr_symbol_mat)))
zero_var_genes <- rowSums(is.na(expr_signature)) > 0
if (any(zero_var_genes)) {
  cat(sum(zero_var_genes), "zero-variance gene(s) dropped before signature construction.\n")
  expr_signature <- expr_signature[!zero_var_genes, ]
}
cat("Shared inference signature matrix:", dim(expr_signature), "\n")

# ------------------------------------------------------------------
# Assemble regulon sources
# ------------------------------------------------------------------
data(dorothea_hs, package = "dorothea")

network_sources <- list(
  dorothea_AC = dorothea_hs %>% dplyr::filter(confidence %in% c("A", "B", "C")),
  dorothea_AB = dorothea_hs %>% dplyr::filter(confidence %in% c("A", "B"))
)

# CollecTRI: fetched from OmnipathR on first run, then CACHED to disk.
# Caching is not a convenience here -- it is a reproducibility
# requirement. A live web resource can change between the run that
# produced the results and the run a reviewer attempts, which would
# silently alter regulon content and therefore every number downstream.
# The cached .rds is the object of record and should be archived with
# the manuscript.
if (file.exists(COLLECTRI_CACHE)) {
  collectri_net <- readRDS(COLLECTRI_CACHE)
  cat("CollecTRI loaded from cache:", COLLECTRI_CACHE, "\n")
} else {
  collectri_net <- tryCatch({
    net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
    saveRDS(net, COLLECTRI_CACHE)
    cat("CollecTRI downloaded and cached (", nrow(net), " edges) on ",
        as.character(Sys.Date()), "\n", sep = "")
    net
  }, error = function(e) {
    warning("CollecTRI could not be retrieved (", conditionMessage(e),
            "). The grid will proceed WITHOUT CollecTRI -- this weakens the ",
            "manuscript's defence against the 'you only tried DoRothEA' ",
            "objection and must be resolved before submission.")
    NULL
  })
}
if (!is.null(collectri_net)) network_sources$collectri <- collectri_net

# Degree-preserving edge-permuted negative controls
if (INCLUDE_EDGEPERM) {
  base_names <- names(network_sources)
  for (i in seq_along(base_names)) {
    nm <- base_names[i]
    network_sources[[paste0(nm, "_edgeperm")]] <-
      permute_network_edges(network_sources[[nm]], seed_offset = 4200L + i)
  }
}

cat("\nRegulon sources assembled:", paste(names(network_sources), collapse = ", "), "\n")

# Normalise schema and filter each network to measured targets
network_grid <- lapply(names(network_sources), function(nm) {
  filt <- filter_network_to_matrix(network_sources[[nm]], expr_signature,
                                   min_targets = MIN_REGULON_SIZE)
  cat(sprintf("  %-24s %5d TFs, %7d edges (after filtering to measured targets)\n",
              nm, nrow(filt$sizes), nrow(filt$network)))
  filt
})
names(network_grid) <- names(network_sources)

regulon_summary <- map_dfr(names(network_grid), function(nm) {
  tibble(regulon = nm,
         n_TFs = nrow(network_grid[[nm]]$sizes),
         n_edges = nrow(network_grid[[nm]]$network),
         median_regulon_size = median(network_grid[[nm]]$sizes$n_targets),
         has_likelihood = any(network_grid[[nm]]$network$likelihood != 1))
})
write_csv(regulon_summary, here::here("results/tables/table_regulon_grid_summary.csv"))

# ------------------------------------------------------------------
# Method runners
#
# Each method is wrapped individually rather than dispatched through
# decoupleR::decouple(), for two reasons: (i) per-method argument
# differences (only viper and wsum consume `likelihood`) can be handled
# explicitly instead of through a nested args list, and (ii) a single
# method failing on a single regulon then degrades to a logged NA for
# that ONE grid cell instead of aborting the entire grid. On a grid
# this size, partial failure must not be fatal.
# ------------------------------------------------------------------
to_matrix <- function(res, score_col = "score") {
  res %>%
    dplyr::select(source, condition, dplyr::all_of(score_col)) %>%
    tidyr::pivot_wider(names_from = condition, values_from = dplyr::all_of(score_col)) %>%
    tibble::column_to_rownames("source") %>%
    as.matrix()
}

run_method <- function(method, network, signature_mat, minsize = MIN_REGULON_SIZE,
                       use_likelihood = TRUE) {

  # `.likelihood` has been deprecated or removed in some decoupleR
  # releases. Rather than let a whole grid cell fail on an argument-name
  # change, the caller retries once with use_likelihood = FALSE; the
  # retry is recorded in the grid log so the manuscript can state which
  # views were likelihood-weighted and which were not, instead of that
  # difference being invisible.
  has_lik <- use_likelihood && any(network$likelihood != 1)

  res <- switch(
    method,

    viper = {
      args <- list(mat = signature_mat, network = network,
                   .source = "source", .target = "target", .mor = "mor",
                   minsize = minsize, verbose = FALSE)
      if (has_lik) args$.likelihood <- "likelihood"
      do.call(decoupleR::run_viper, args)
    },

    ulm = decoupleR::run_ulm(
      mat = signature_mat, network = network,
      .source = "source", .target = "target", .mor = "mor", minsize = minsize
    ),

    mlm = decoupleR::run_mlm(
      mat = signature_mat, network = network,
      .source = "source", .target = "target", .mor = "mor", minsize = minsize
    ),

    norm_wsum = {
      if (WSUM_PERMUTATIONS <= 0) stop("norm_wsum disabled by config")
      set_pipeline_seed(offset = 5100)
      args <- list(mat = signature_mat, network = network,
                   .source = "source", .target = "target", .mor = "mor",
                   times = WSUM_PERMUTATIONS, seed = PIPELINE_SEED,
                   minsize = minsize)
      if (has_lik) args$.likelihood <- "likelihood"
      do.call(decoupleR::run_wsum, args)
    },

    stop("Unknown method: ", method)
  )

  # decoupleR returns a `statistic` column when a call emits several
  # related statistics (wsum / corr_wsum / norm_wsum). Select the one
  # actually requested; do not average or silently take the first.
  if ("statistic" %in% colnames(res) && dplyr::n_distinct(res$statistic) > 1) {
    wanted <- if (method == "norm_wsum") "norm_wsum" else method
    if (!wanted %in% unique(res$statistic)) {
      stop("Requested statistic '", wanted, "' not present in decoupleR output. ",
           "Available: ", paste(unique(res$statistic), collapse = ", "))
    }
    res <- dplyr::filter(res, statistic == wanted)
  }

  score_col <- intersect(c("score", "statistic_value", "estimate"), colnames(res))[1]
  if (is.na(score_col)) {
    stop("No recognised score column in decoupleR output for method '", method,
         "'. Columns present: ", paste(colnames(res), collapse = ", "))
  }
  to_matrix(res, score_col)
}

# ------------------------------------------------------------------
# Run the grid
# ------------------------------------------------------------------
cat("\n================ Running TF activity grid ================\n")

# ORDERING: METHOD-MAJOR, in priority order. expand_grid() varies its
# LAST argument fastest, so listing method first means the grid runs
# viper across every regulon, then ulm across every regulon, and only
# then the more speculative methods.
#
# This matters because the grid can run for hours and may have to be
# stopped. Regulon-major ordering would complete DoRothEA in all four
# methods before touching CollecTRI at all -- so an interrupted run
# would be missing the single most important regulon source. Under
# method-major ordering an interrupted run still has a complete,
# defensible {all regulons} x {viper, ulm} grid, which is the part the
# manuscript actually rests on.
METHOD_PRIORITY <- c("viper", "ulm", "mlm", "norm_wsum")
methods_ordered <- c(intersect(METHOD_PRIORITY, METHODS_REQUESTED),
                     setdiff(METHODS_REQUESTED, METHOD_PRIORITY))

grid_spec <- tidyr::expand_grid(
  method  = methods_ordered,
  regulon = names(network_grid)
)

cat("Grid:", nrow(grid_spec), "cells --", length(names(network_grid)),
    "regulons x", length(methods_ordered), "methods\n")
cat("Method order (interrupt-safe):", paste(methods_ordered, collapse = " -> "), "\n\n")

# Resume support: a completed view is checkpointed immediately, so a
# crash or a deliberate stop does not throw away hours of VIPER runs.
GRID_CHECKPOINT <- here::here("data/processed/tf_activity_grid_partial.rds")
activity_grid <- if (file.exists(GRID_CHECKPOINT)) {
  ck <- readRDS(GRID_CHECKPOINT)
  cat("Resuming from checkpoint:", length(ck), "view(s) already complete --",
      paste(names(ck), collapse = ", "), "\n")
  cat("Delete", basename(GRID_CHECKPOINT), "to force a clean re-run.\n\n")
  ck
} else list()
grid_log <- list()

for (i in seq_len(nrow(grid_spec))) {
  reg_nm <- grid_spec$regulon[i]
  mth    <- grid_spec$method[i]
  view_id <- paste(reg_nm, mth, sep = "__")

  if (view_id %in% names(activity_grid)) {
    cat(sprintf("  [skip]   %-34s already in checkpoint\n", view_id))
    grid_log[[view_id]] <- tibble(
      view_id = view_id, regulon = reg_nm, method = mth, status = "ok",
      n_TFs = nrow(activity_grid[[view_id]]), runtime_min = NA_real_,
      message = "restored from checkpoint")
    next
  }

  t0 <- Sys.time()
  attempt <- function(use_lik) {
    m <- run_method(mth, network_grid[[reg_nm]]$network, expr_signature,
                    use_likelihood = use_lik)
    # Align to the canonical sample order used everywhere downstream.
    m <- m[, rownames(sample_metadata), drop = FALSE]
    stopifnot(identical(colnames(m), rownames(sample_metadata)))
    m
  }
  out <- tryCatch({
    list(matrix = attempt(TRUE), status = "ok", message = NA_character_)
  }, error = function(e1) {
    tryCatch({
      list(matrix = attempt(FALSE), status = "ok",
           message = paste0("retried without .likelihood after: ", conditionMessage(e1)))
    }, error = function(e2) {
      list(matrix = NULL, status = "failed", message = conditionMessage(e2))
    })
  })
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  if (out$status == "ok") {
    activity_grid[[view_id]] <- out$matrix
    saveRDS(activity_grid, GRID_CHECKPOINT)     # checkpoint after every view
    cat(sprintf("  [ok]     %-34s %4d TFs  (%.1f min)\n",
                view_id, nrow(out$matrix), elapsed))
  } else {
    cat(sprintf("  [FAILED] %-34s %s\n", view_id, out$message))
  }

  grid_log[[view_id]] <- tibble(
    view_id = view_id, regulon = reg_nm, method = mth,
    status = out$status,
    n_TFs = if (out$status == "ok") nrow(out$matrix) else NA_integer_,
    runtime_min = elapsed,
    message = out$message
  )
}

# ------------------------------------------------------------------
# Consensus view: rank-based consensus across the successful methods
# within each NON-permuted regulon.
#
# Computed on per-TF rank-normalised scores, because the methods live
# on incomparable scales (VIPER normalised enrichment vs a t-like ulm
# statistic vs a permutation z-score). Averaging raw scores across
# methods would silently let whichever method has the widest numeric
# range dominate the consensus.
# ------------------------------------------------------------------
rank_normalise <- function(m) {
  t(apply(m, 1, function(r) {
    if (all(is.na(r))) return(r)
    stats::qnorm((rank(r, na.last = "keep") - 0.5) / sum(!is.na(r)))
  }))
}

for (reg_nm in names(network_grid)) {
  members <- grep(paste0("^", reg_nm, "__"), names(activity_grid), value = TRUE)
  if (length(members) < 2) next
  common_tfs <- Reduce(intersect, lapply(activity_grid[members], rownames))
  if (length(common_tfs) < 10) next
  stacked <- lapply(members, function(v) rank_normalise(activity_grid[[v]][common_tfs, , drop = FALSE]))
  consensus_mat <- Reduce(`+`, stacked) / length(stacked)
  view_id <- paste(reg_nm, "consensus", sep = "__")
  activity_grid[[view_id]] <- consensus_mat
  grid_log[[view_id]] <- tibble(
    view_id = view_id, regulon = reg_nm, method = "consensus",
    status = "ok", n_TFs = nrow(consensus_mat), runtime_min = 0,
    message = paste("rank-consensus of:", paste(members, collapse = "; "))
  )
  cat(sprintf("  [ok]     %-34s %4d TFs  (consensus of %d methods)\n",
              view_id, nrow(consensus_mat), length(members)))
}

grid_log_table <- bind_rows(grid_log)
write_csv(grid_log_table, here::here("results/tables/table_tf_activity_grid_log.csv"))

cat("\nGrid complete:", sum(grid_log_table$status == "ok"), "of",
    nrow(grid_log_table), "views produced successfully.\n")

if (sum(grid_log_table$status == "ok") < 4) {
  stop("Fewer than 4 TF-activity views were produced. The grid is too sparse ",
       "to support a defensible multi-configuration claim -- resolve the ",
       "failures above before proceeding to script 07.")
}

# ------------------------------------------------------------------
# Cross-view concordance
#
# Interpretation matters here, and it cuts both ways:
#   HIGH concordance across methods => the downstream null result is a
#     property of TF activity as a construct, not of one scoring
#     algorithm. That STRENGTHENS a negative conclusion.
#   LOW concordance => the methods disagree about what TF activity even
#     is, so "TF activity adds nothing" would be an overreach and the
#     claim must be stated per-method.
# Either way this table must be reported, because it determines how
# broadly the paper is entitled to generalise.
# ------------------------------------------------------------------
ok_views <- names(activity_grid)
concordance_pairs <- combn(ok_views, 2, simplify = FALSE)

view_concordance <- map_dfr(concordance_pairs, function(pr) {
  a <- activity_grid[[pr[1]]]; b <- activity_grid[[pr[2]]]
  shared <- intersect(rownames(a), rownames(b))
  if (length(shared) < 10) {
    return(tibble(view_a = pr[1], view_b = pr[2], n_shared_TFs = length(shared),
                  median_per_TF_spearman = NA_real_))
  }
  rhos <- vapply(shared, function(tfx) {
    suppressWarnings(cor(a[tfx, ], b[tfx, ], method = "spearman"))
  }, numeric(1))
  tibble(view_a = pr[1], view_b = pr[2], n_shared_TFs = length(shared),
         median_per_TF_spearman = median(rhos, na.rm = TRUE))
})
write_csv(view_concordance, here::here("results/tables/table_tf_view_concordance.csv"))

cat("\nCross-view concordance (median per-TF Spearman), most similar pairs:\n")
print(view_concordance %>% arrange(desc(median_per_TF_spearman)) %>% head(10))
cat("\nLeast similar pairs:\n")
print(view_concordance %>% arrange(median_per_TF_spearman) %>% head(10))

# Real-vs-edge-permuted concordance is the sharpest single number in
# this script. If a real regulon's activity matrix correlates highly
# with its own degree-matched edge-permuted control, then the "TF
# activity" signal is dominated by regulon SIZE and gene-set averaging
# rather than by curated regulatory content -- which would explain a
# downstream null far more informatively than "integration did not
# help", and would itself be a publishable methodological observation.
edgeperm_check <- view_concordance %>%
  filter(str_detect(view_b, "_edgeperm__") |  str_detect(view_a, "_edgeperm__")) %>%
  mutate(
    base_a = str_remove(view_a, "_edgeperm"),
    base_b = str_remove(view_b, "_edgeperm")
  ) %>%
  filter(base_a == base_b) %>%
  dplyr::select(view_a, view_b, n_shared_TFs, median_per_TF_spearman)

if (nrow(edgeperm_check) > 0) {
  cat("\n=== Real regulon vs its own degree-matched edge-permuted control ===\n")
  print(edgeperm_check)
  write_csv(edgeperm_check, here::here("results/tables/table_real_vs_edgepermuted_concordance.csv"))
}

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------
saveRDS(activity_grid, here::here("data/processed/tf_activity_grid.rds"))
# The grid completed, so the partial checkpoint is now redundant.
if (file.exists(GRID_CHECKPOINT)) unlink(GRID_CHECKPOINT)

# BACKWARD COMPATIBILITY: scripts 06 and 12 read these exact filenames.
# The primary view is unchanged (DoRothEA A-C + VIPER), so every
# previously reported downstream number remains reproducible bit-for-bit
# and this revision is strictly additive rather than a replacement.
primary_view <- "dorothea_AC__viper"
if (primary_view %in% names(activity_grid)) {
  saveRDS(activity_grid[[primary_view]],
          here::here("data/processed/tf_activity_viper_AC_primary.rds"))
  cat("\nBackward-compatible primary view written:", primary_view, "\n")
} else {
  warning("Primary view '", primary_view, "' failed. Downstream scripts 06/12 ",
          "expect data/processed/tf_activity_viper_AC_primary.rds and will not run.")
}
if ("dorothea_AB__viper" %in% names(activity_grid)) {
  saveRDS(activity_grid[["dorothea_AB__viper"]],
          here::here("data/processed/tf_activity_viper_AB_sensitivity.rds"))
}
saveRDS(network_grid[["dorothea_AC"]]$sizes, here::here("data/processed/tf_size_AC.rds"))
write_csv(network_grid[["dorothea_AC"]]$sizes,
          here::here("results/tables/table_tf_regulon_sizes_AC.csv"))

# ------------------------------------------------------------------
# Heatmap of the primary view (display only)
# ------------------------------------------------------------------
if (primary_view %in% names(activity_grid)) {
  tf_primary <- activity_grid[[primary_view]]
  tf_var    <- apply(tf_primary, 1, var)
  top50_tfs <- names(sort(tf_var, decreasing = TRUE))[1:min(50, length(tf_var))]
  annotation_col <- sample_metadata %>%
    dplyr::select(PAM50) %>%
    mutate(PAM50 = ifelse(is.na(PAM50) | PAM50 == "", "Unknown", PAM50))

  png(here::here("results/figures/figure_tf_activity_top50_heatmap.png"),
      width = 1600, height = 1200, res = 150)
  pheatmap(
    tf_primary[top50_tfs, ], scale = "row",
    annotation_col = annotation_col, show_colnames = FALSE,
    main = "Top-50-variance TF activity (VIPER, DoRothEA A-C), z-scored per TF for display only"
  )
  dev.off()
}

# The design note is generated FROM THE ACTUAL RUN, not written as a
# fixed block of prose. A hard-coded note describing four methods and a
# CollecTRI arm would silently become false the moment the CONFIG block
# is trimmed -- and a methods description that does not match what the
# code executed is exactly the kind of discrepancy that destroys trust
# in a thesis when an examiner checks it against the logs.
methods_run  <- sort(unique(grid_log_table$method[grid_log_table$status == "ok"]))
regulons_run <- sort(unique(grid_log_table$regulon[grid_log_table$status == "ok"]))
edgeperm_run <- any(grepl("_edgeperm$", regulons_run))
collectri_run <- any(grepl("^collectri", regulons_run))

writeLines(c(
  "DESIGN NOTE: TF-activity inference is a GRID, not a single view.",
  "======================================================================",
  paste0("Generated from the actual run on ", as.character(Sys.Date()), "."),
  "",
  "Rationale: a null added-value result obtained from ONE regulon and ONE",
  "inference method cannot be distinguished from a misconfiguration. This",
  "script therefore produces a grid of regulon x method views, and script",
  "07 evaluates added value across them.",
  "",
  paste0("REGULON SOURCES ACTUALLY USED (", length(regulons_run), "): ",
         paste(regulons_run, collapse = ", ")),
  paste0("INFERENCE METHODS ACTUALLY USED (", length(methods_run), "): ",
         paste(methods_run, collapse = ", ")),
  paste0("TOTAL VIEWS PRODUCED: ", length(activity_grid)),
  "",
  if (edgeperm_run)
    paste0("EDGE-PERMUTED CONTROLS: INCLUDED. These keep every TF's out-degree ",
           "and the mor/likelihood distribution but randomise WHICH genes each ",
           "TF regulates, separating 'curated regulatory biology contributes' ",
           "from 'averaging arbitrary gene sets of the same size contributes'. ",
           "This is a strictly stronger control than the sample-permutation ",
           "control in script 07; the two are NOT interchangeable.")
  else
    paste0("EDGE-PERMUTED CONTROLS: *** NOT INCLUDED *** (INCLUDE_EDGEPERM was ",
           "FALSE). Without them the analysis cannot distinguish a benefit ",
           "arising from curated regulatory content from one arising merely ",
           "from averaging gene sets of the same size. Script 07's default ",
           "VIEWS_FOR_BOOTSTRAP also names edge-permuted views, which will be ",
           "silently dropped. THIS MUST BE RESOLVED BEFORE SUBMISSION."),
  "",
  if (collectri_run)
    "REGULON SOURCE VARIED: yes -- both DoRothEA and CollecTRI are present."
  else
    paste0("REGULON SOURCE: *** DoRothEA ONLY *** -- CollecTRI is absent. The ",
           "grid therefore varies inference method and confidence tier but NOT ",
           "regulon source, and does not answer the 'you only tried DoRothEA' ",
           "objection. State this as an explicit limitation."),
  "",
  "All methods consume the IDENTICAL per-gene z-scored signature, so grid",
  "cells differ only in the factor being varied.",
  "",
  "The primary view (dorothea_AC__viper) is unchanged from the previous",
  "revision and is still written to its original filename, so all",
  "previously reported downstream results remain exactly reproducible.",
  "",
  "See also: table_tf_activity_grid_log.csv (what ran, how long, what",
  "failed), table_tf_view_concordance.csv (how much the views agree --",
  "this determines how broadly any conclusion may be generalised) and",
  "table_real_vs_edgepermuted_concordance.csv (whether the signal is",
  "regulatory content or gene-set averaging)."
), here::here("results/tables/NOTE_tf_activity_grid_design.txt"))

if (!edgeperm_run || !collectri_run) {
  warning("Grid is INCOMPLETE: ",
          if (!edgeperm_run) "no edge-permuted controls. " else "",
          if (!collectri_run) "no CollecTRI. " else "",
          "See NOTE_tf_activity_grid_design.txt. Script 07 will bootstrap ",
          "fewer views than its pre-specified list.")
}

log_session_info(script_name, key_packages = c("dorothea", "decoupleR", "OmnipathR"))
cat("\n✓ 05_tf_activity.R complete.",
    length(activity_grid), "TF-activity views written to tf_activity_grid.rds\n")

close_logger(log_con, script_name)
