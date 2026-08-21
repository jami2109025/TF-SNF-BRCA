# ============================================================
# R/utils.R  (publication-ready revision)
# Shared utilities for the BRCA SNF-TF subtyping pipeline.
#
# This revision fixes the following issues identified in the pipeline
# audit (kept as a durable changelog so provenance of each helper is
# traceable for peer review):
#   - map_to_symbol_dedup() now accepts an explicit `variance_reference`
#     matrix so the SAME underlying values (always the VST matrix) are
#     used to break ties between duplicate Ensembl IDs, regardless of
#     which matrix (VST or raw counts) is ultimately being deduplicated.
#     Previously, script 05 tie-broke on VST variance and script 11 on
#     raw-count variance, which could silently select a DIFFERENT
#     representative Ensembl ID for the same gene symbol in each script.
#   - pick_best_covariate() added: chooses a clinical covariate column
#     by DATA COMPLETENESS, not by first-match-in-a-candidate-list
#     order, so an external cohort's Cox model is not silently
#     underpowered by picking a mostly-missing column.
#   - cohens_d_one_vs_rest() added: a genuine pooled-SD, one-vs-rest
#     standardized effect size (matching the InCluster/Other framing
#     used consistently elsewhere in the pipeline, e.g. DESeq2 in
#     script 11), replacing the previous cluster-vs-grand-mean /
#     within-cluster-SD statistic that was labeled "Cohen's-d-like"
#     but did not reduce to the standard two-group formula.
#   - test_cluster_batch_association() added: a formal, reusable test
#     of cluster assignment against a technical batch proxy (plate,
#     sequencing center), closing the audit's "batch proxies computed
#     but never tested against final clusters" gap.
#   - require_columns(), align_samples(), cluster_factor(),
#     resolve_clinical_columns(), logging/seeding infrastructure are
#     carried over unchanged from the prior revision; they were
#     verified correct in the audit.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

# ------------------------------------------------------------------
# 1. Reproducibility / logging infrastructure
# ------------------------------------------------------------------

PIPELINE_SEED <- 20240601L

#' Set the pipeline-wide seed, optionally offset for a specific
#' stochastic sub-step. Every script/sub-step calls this explicitly
#' instead of an ad hoc set.seed(), so the seed provenance for any
#' given random draw is always traceable to one line of code.
set_pipeline_seed <- function(offset = 0L) {
  set.seed(PIPELINE_SEED + offset)
}

#' Lightweight script-level logger: console + per-script log file.
init_logger <- function(script_name) {
  dir.create(here::here("logs"), showWarnings = FALSE, recursive = TRUE)
  log_path <- here::here("logs", paste0(script_name, ".log"))
  con <- file(log_path, open = "wt")
  sink(con, split = TRUE)
  sink(con, type = "message", append = TRUE)
  cat(sprintf("==== %s started: %s ====\n", script_name, Sys.time()))
  cat("Pipeline seed:", PIPELINE_SEED, "\n")
  invisible(con)
}

close_logger <- function(con, script_name) {
  cat(sprintf("==== %s finished: %s ====\n", script_name, Sys.time()))
  sink(type = "message")
  sink()
  close(con)
}

#' Persist sessionInfo() + a tidy key-package-version table after every
#' script, so a package upgrade between two pipeline runs is detectable
#' after the fact rather than only inferred from a numerical mismatch.
log_session_info <- function(script_name, key_packages = NULL) {
  dir.create(here::here("logs", "sessioninfo"), showWarnings = FALSE, recursive = TRUE)
  writeLines(
    capture.output(sessionInfo()),
    here::here("logs", "sessioninfo", paste0(script_name, "_sessionInfo.txt"))
  )
  if (!is.null(key_packages)) {
    versions <- tibble(
      package = key_packages,
      version = vapply(key_packages, function(p) {
        tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
      }, character(1)),
      script  = script_name,
      date    = as.character(Sys.Date())
    )
    write_csv(
      versions,
      here::here("logs", "sessioninfo", paste0(script_name, "_package_versions.csv"))
    )
  }
}

#' Fail loudly and informatively instead of silently proceeding with NA.
require_columns <- function(df, cols, context = "") {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop(
      "[", context, "] Missing required column(s): ",
      paste(missing, collapse = ", "),
      ". Available columns: ", paste(colnames(df), collapse = ", ")
    )
  }
  invisible(TRUE)
}

#' Idempotent directory scaffold, callable from any script in isolation.
ensure_dirs <- function() {
  folders <- c(
    "data/raw", "data/processed", "data/external",
    "results/tables", "results/figures", "results/objects",
    "logs", "logs/sessioninfo"
  )
  for (f in folders) dir.create(here::here(f), recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------------
# 2. Single source of truth for clinical column names
# ------------------------------------------------------------------

CLINICAL_COLUMN_MAP <- list(
  patient_id             = "patient_id",
  vital_status           = "paper_vital_status",
  days_to_death          = "paper_days_to_death",
  days_to_last_followup  = "paper_days_to_last_followup",
  age_years              = "paper_age_at_initial_pathologic_diagnosis",
  stage                  = "paper_pathologic_stage",
  pathology              = "paper_BRCA_Pathology",
  grade                  = "paper_Tumor_Grade",
  er_status              = "paper_ER_Status",
  pr_status              = "paper_PR_Status",
  her2_status            = "paper_HER2_Final_Status"
)

#' Resolve and standardize the authoritative clinical columns onto a
#' metadata data.frame. Called once (script 03); every script after
#' that consumes ONLY these resolved names.
resolve_clinical_columns <- function(df) {

  resolved <- tibble(.rows = nrow(df))

  for (concept in names(CLINICAL_COLUMN_MAP)) {
    src_col <- CLINICAL_COLUMN_MAP[[concept]]
    if (src_col %in% colnames(df)) {
      resolved[[concept]] <- df[[src_col]]
    } else {
      warning(
        "resolve_clinical_columns(): expected source column '", src_col,
        "' for concept '", concept, "' not found; filling with NA. ",
        "Check CLINICAL_COLUMN_MAP in R/utils.R against the columns ",
        "actually returned by GDCquery_clinic()/colData() for this ",
        "TCGAbiolinks version."
      )
      resolved[[concept]] <- NA
    }
  }

  resolved <- resolved %>%
    mutate(
      stage_simple = case_when(
        is.na(stage) ~ NA_character_,
        grepl("^Stage_?IV",  stage, ignore.case = TRUE) ~ "Stage IV",
        grepl("^Stage_?III", stage, ignore.case = TRUE) ~ "Stage III",
        grepl("^Stage_?II",  stage, ignore.case = TRUE) ~ "Stage II",
        grepl("^Stage_?I\\b|^Stage_?I$", stage, ignore.case = TRUE) ~ "Stage I",
        TRUE ~ NA_character_
      ),
      stage_simple = factor(stage_simple,
                             levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
      pathology_clean = ifelse(
        pathology %in% c("IDC", "ILC", "Mixed", "Other"),
        pathology, NA_character_
      ),
      age_years = suppressWarnings(as.numeric(age_years))
      # days_to_death = suppressWarnings(as.numeric(days_to_death)),
      # days_to_last_followup = suppressWarnings(as.numeric(days_to_last_followup))
    )

  resolved
}

#' Data-completeness-based covariate selection for an EXTERNAL cohort
#' whose clinical schema doesn't match TCGA's curated fields (e.g.
#' METABRIC). Fixes the audit finding that first-match-in-a-list
#' selection could silently choose a mostly-missing column over a more
#' complete alternative, underpowering downstream Cox models.
#'
#' @param candidates character vector of candidate column names, in
#'   NO particular priority order (completeness decides, not order).
#' @param data data.frame to search.
#' @param min_complete_frac warn if the best available column is still
#'   below this completeness fraction.
pick_best_covariate <- function(candidates, data, min_complete_frac = 0.5) {
  candidates <- candidates[candidates %in% colnames(data)]
  if (length(candidates) == 0) return(NA_character_)
  completeness <- vapply(candidates, function(cn) {
    v <- data[[cn]]
    mean(!is.na(v) & !(is.character(v) & v == ""))
  }, numeric(1))
  best <- names(completeness)[which.max(completeness)]
  if (completeness[best] < min_complete_frac) {
    warning(
      "pick_best_covariate(): best available column '", best, "' is only ",
      round(100 * completeness[best], 1), "% complete (candidates checked: ",
      paste(candidates, collapse = ", "), "). Proceeding, but this may ",
      "materially reduce complete-case N in downstream models."
    )
  }
  best
}

# ------------------------------------------------------------------
# 3. Cluster-label handling
# ------------------------------------------------------------------

#' Dynamic-level cluster factor -- NEVER hardcode levels, since k varies
#' (k=2 and k=3 are both run throughout this pipeline).
cluster_factor <- function(snf_cluster_vector) {
  lv <- sort(unique(stats::na.omit(snf_cluster_vector)))
  factor(snf_cluster_vector, levels = lv)
}

#' Formal test of SNF cluster assignment against a technical batch
#' proxy variable (e.g. TCGA barcode-derived plate/sequencing center).
#' Closes the audit finding that batch proxies were computed and
#' plotted (script 04) but never formally tested against final cluster
#' labels once those existed (script 06+).
#'
#' Uses chi-square (or Monte-Carlo-simulated Fisher's exact when
#' expected cell counts are small) on the cluster x batch-proxy
#' contingency table, exactly mirroring the test-selection logic used
#' for PAM50/clinical concordance elsewhere in this pipeline, for
#' methodological consistency.
test_cluster_batch_association <- function(cluster_vector, batch_vector, seed_offset = 0L) {
  tab <- table(cluster_vector, batch_vector, useNA = "no")
  tab <- tab[, colSums(tab) > 0, drop = FALSE]
  tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(list(test = "insufficient data", statistic = NA_real_, p_value = NA_real_, table = tab))
  }
  chi <- suppressWarnings(chisq.test(tab))
  if (any(chi$expected < 5)) {
    set_pipeline_seed(offset = seed_offset)
    ft <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
    list(test = "Fisher exact (simulated)", statistic = NA_real_, p_value = ft$p.value, table = tab)
  } else {
    list(test = "Chi-square", statistic = unname(chi$statistic), p_value = chi$p.value, table = tab)
  }
}

# ------------------------------------------------------------------
# 4. Sample-matching helper
# ------------------------------------------------------------------

#' Intersect sample IDs across an arbitrary number of NAMED objects
#' (matrices by colname, data.frames by rownames), return re-ordered
#' identically, with a hard stopifnot() verification.
align_samples <- function(...) {
  objs <- list(...)
  obj_names <- names(objs)
  if (is.null(obj_names) || any(obj_names == "")) {
    stop("align_samples(): every argument must be named, e.g. ",
         "align_samples(expr = mat1, tf = mat2, meta = df1)")
  }

  get_ids <- function(x) if (is.matrix(x)) colnames(x) else rownames(x)

  id_list <- lapply(objs, get_ids)
  common  <- Reduce(intersect, id_list)
  if (length(common) == 0) {
    stop("align_samples(): no common sample IDs found across: ",
         paste(obj_names, collapse = ", "))
  }

  aligned <- lapply(objs, function(x) {
    if (is.matrix(x)) x[, common, drop = FALSE] else x[common, , drop = FALSE]
  })
  names(aligned) <- obj_names

  ref <- get_ids(aligned[[1]])
  for (nm in obj_names[-1]) {
    stopifnot(identical(get_ids(aligned[[nm]]), ref))
  }

  cat("align_samples(): ", length(common), " common samples across {",
      paste(obj_names, collapse = ", "), "}\n", sep = "")

  aligned
}

# ------------------------------------------------------------------
# 5. Gene ID <-> symbol mapping, deduplication
# ------------------------------------------------------------------

#' Map Ensembl IDs (rownames of `mat`) to gene SYMBOL, collapsing
#' duplicate symbols by keeping the row with the highest variance IN
#' A CALLER-SUPPLIED REFERENCE MATRIX (`variance_reference`), which
#' defaults to `mat` itself but should be passed explicitly as the
#' SAME matrix (the VST matrix) everywhere in the pipeline that this
#' function is called, so the tie-break rule is scale-consistent and
#' the SAME representative Ensembl ID is chosen for a given duplicate
#' symbol regardless of what matrix (VST or raw counts) is ultimately
#' being deduplicated and returned.
#'
#' Fixes the audit finding that script 05 (tie-breaking on VST
#' variance) and script 11 (tie-breaking on raw-count variance, via
#' the same function but a different input matrix) could silently
#' select DIFFERENT Ensembl IDs for the same gene symbol.
#'
#' @param mat matrix to deduplicate and return (rownames = Ensembl IDs,
#'   possibly versioned e.g. "ENSG00000141510.16").
#' @param variance_reference matrix used ONLY to compute the
#'   tie-break variance; must share the same Ensembl IDs (rownames) as
#'   `mat` after version-stripping. Pass the VST matrix explicitly.
#' @param org_db an AnnotationDbi OrgDb object (e.g. org.Hs.eg.db).
map_to_symbol_dedup <- function(mat, variance_reference = mat, org_db = org.Hs.eg.db::org.Hs.eg.db) {

  ensembl_clean     <- sub("\\..*$", "", rownames(mat))
  ensembl_clean_ref <- sub("\\..*$", "", rownames(variance_reference))

  symbol_map <- AnnotationDbi::select(
    org_db,
    keys    = unique(ensembl_clean),
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "ENTREZID")
  ) %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::distinct(ENSEMBL, SYMBOL, ENTREZID, .keep_all = FALSE)

  # Variance is always computed from variance_reference, keyed by the
  # SAME (version-stripped) Ensembl ID, then joined onto `mat`'s rows.
  ref_variance <- tibble(
    ENSEMBL = ensembl_clean_ref,
    gene_variance = apply(variance_reference, 1, var)
  ) %>%
    dplyr::group_by(ENSEMBL) %>%
    dplyr::summarise(gene_variance = max(gene_variance, na.rm = TRUE), .groups = "drop")

  df <- as.data.frame(mat)
  df$ENSEMBL <- ensembl_clean

  df_symbol <- df %>%
    dplyr::inner_join(symbol_map, by = "ENSEMBL", relationship = "many-to-many") %>%
    dplyr::left_join(ref_variance, by = "ENSEMBL") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(order_by = gene_variance, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  out_mat <- df_symbol %>%
    dplyr::select(dplyr::all_of(colnames(mat))) %>%
    as.matrix()
  rownames(out_mat) <- df_symbol$SYMBOL

  list(
    matrix = out_mat,
    map    = df_symbol %>% dplyr::select(SYMBOL, ENSEMBL, ENTREZID)
  )
}

# ------------------------------------------------------------------
# 6. Standardized effect sizes
# ------------------------------------------------------------------

#' Genuine one-vs-rest Cohen's d with pooled SD, mirroring the
#' InCluster/Other framing used consistently elsewhere in this
#' pipeline (e.g. DESeq2 contrasts in script 11). Replaces the
#' previous "cluster mean minus grand mean, divided by within-cluster
#' SD" statistic, which was labeled Cohen's-d-like but does not reduce
#' to the standard two-group formula for k > 2, and could rank TFs
#' differently -- especially with unequal cluster sizes.
#'
#' @param values numeric vector (e.g. TF activity for one TF, all samples)
#' @param group_vector cluster/group assignment, same length as values
#' @param target_group the group treated as "in-cluster"
cohens_d_one_vs_rest <- function(values, group_vector, target_group) {
  in_grp  <- values[group_vector == target_group]
  out_grp <- values[group_vector != target_group]
  n1 <- sum(!is.na(in_grp)); n2 <- sum(!is.na(out_grp))
  if (n1 < 2 || n2 < 2) return(NA_real_)
  m1 <- mean(in_grp, na.rm = TRUE);  s1 <- sd(in_grp, na.rm = TRUE)
  m2 <- mean(out_grp, na.rm = TRUE); s2 <- sd(out_grp, na.rm = TRUE)
  pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(pooled_sd) || pooled_sd <= 0) return(NA_real_)
  (m1 - m2) / pooled_sd
}

cat("R/utils.R loaded (publication-ready revision). Pipeline seed =", PIPELINE_SEED, "\n")
