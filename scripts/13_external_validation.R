# ============================================================
# 13_external_validation.R  (publication-ready revision, FINAL+PH-FIX)
# External validation: project METABRIC onto TCGA SNF centroids

suppressPackageStartupMessages({
  library(here)
  library(org.Hs.eg.db)
  library(survival)
  library(survminer)
  library(ggpubr)
  library(survRM2)      # NEW: PH-assumption-free RMST comparisons
  library(tidyverse)
})

source(here::here("R", "utils.R"))

script_name <- "13_external_validation"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

metabric_expr_path <- here::here("data/external/metabric_expression.txt")
metabric_clin_path <- here::here("data/external/metabric_clinical.txt")
if (!file.exists(metabric_expr_path) || !file.exists(metabric_clin_path)) {
  stop(
    "METABRIC files not found at data/external/metabric_expression.txt and ",
    "data/external/metabric_clinical.txt. These must be obtained manually ",
    "(e.g. from cBioPortal) and placed there before running this script."
  )
}

vst_top5000           <- readRDS(here::here("data/processed/vst_top5000_genes.rds"))
vst_matrix_all_genes  <- readRDS(here::here("data/processed/vst_matrix_all_genes.rds"))
final_labels_obj      <- readRDS(here::here("results/objects/final_snf_cluster_labels.rds"))

cross_platform_note <- c(
  "METHODOLOGICAL CAVEAT: cross-platform / cross-population comparability",
  "==========================================================================",
  "TCGA expression is RNA-seq; METABRIC expression is typically microarray",
  "(Illumina HT-12). Both cohorts are z-scored INDEPENDENTLY (within-cohort,",
  "per-gene) before nearest-centroid correlation classification, which is",
  "standard practice and avoids cross-platform batch-leakage that joint",
  "normalization would introduce. Independent z-scoring corrects only for",
  "per-gene SCALE and LOCATION differences -- it does NOT correct for",
  "platform-specific non-linear measurement differences, nor for population",
  "differences between TCGA (US-based, more ethnically diverse) and METABRIC",
  "(UK-based) cohorts. Classification confidence and replication strength",
  "should be interpreted with this caveat in mind.",
  "",
  "COVARIATE-MATCHING NOTE (added in this revision): the TCGA discovery-cohort",
  "multivariable Cox model adjusts for age AND stage. This METABRIC clinical",
  "export (cBioPortal brca_metabric clinical-patient file) does NOT contain a",
  "tumor stage or grade column at all -- the METABRIC multivariable model here",
  "is therefore adjusted for age and ER (IHC) status ONLY. This is a genuine",
  "gap in covariate availability between the two cohorts, not a matched",
  "replication of the discovery-cohort adjustment set, and should be",
  "interpreted as such.",
  "",
  "PROPORTIONAL-HAZARDS CAVEAT (added in this revision): cox.zph() diagnostics",
  "show significant non-proportional hazards for the assigned_SNF_cluster term",
  "(and covariates) in both the k=2 and k=3 multivariable Cox models -- visibly",
  "consistent with the crossing Kaplan-Meier curves in both external-validation",
  "KM figures. The single hazard ratios reported in",
  "table_metabric_cox_multivariable_k2_k3.csv are therefore time-averaged",
  "effects, not constant relative risks. Schoenfeld-residual plots, a",
  "time-split Cox model, and RMST pairwise comparisons are provided as",
  "PH-assumption-free (or PH-diagnostic) sensitivity analyses; see",
  "figure_metabric_schoenfeld_k*.png, table_metabric_cox_timesplit_k2_k3.csv,",
  "and table_metabric_rmst_k2_k3.csv."
)
writeLines(cross_platform_note, here::here("results/tables/NOTE_cross_platform_caveat.txt"))
cat(cross_platform_note, sep = "\n"); cat("\n")

# -----------------------------
# Map TCGA top-5000 feature set from Ensembl IDs to gene symbols, using
# the SAME variance-reference convention as scripts 05/11, before doing
# anything cross-cohort. Comparing Ensembl-ID rownames against
# METABRIC's gene-symbol rownames would otherwise silently intersect to
# zero genes.
# -----------------------------
mapped_5000 <- map_to_symbol_dedup(
  vst_top5000,
  variance_reference = vst_matrix_all_genes,
  org_db = org.Hs.eg.db
)
vst_top5000_symbol <- mapped_5000$matrix
cat("vst_top5000 mapped from", nrow(vst_top5000), "Ensembl IDs to",
    nrow(vst_top5000_symbol), "unique gene symbols\n")

# -----------------------------
# Build TCGA per-cluster centroids (z-scored within TCGA, by gene)
# -----------------------------
tcga_z <- t(scale(t(vst_top5000_symbol)))

build_centroids <- function(cluster_labels_k, expr_z) {
  expr_z_aligned <- expr_z[, cluster_labels_k$full_barcode]
  stopifnot(identical(colnames(expr_z_aligned), cluster_labels_k$full_barcode))
  sapply(split(cluster_labels_k$full_barcode, cluster_labels_k$SNF_cluster), function(ids) {
    rowMeans(expr_z_aligned[, ids, drop = FALSE])
  })
}

centroids_k2 <- build_centroids(final_labels_obj$k2, tcga_z)
centroids_k3 <- build_centroids(final_labels_obj$k3, tcga_z)
cat("TCGA centroids built: k2 =", ncol(centroids_k2), "clusters, k3 =", ncol(centroids_k3), "clusters\n")

# -----------------------------
# Load + clean METABRIC expression, with an explicit safety check that
# every remaining numeric column is actually a sample (any leftover
# numeric metadata column would otherwise be silently averaged in as if
# it were a sample's expression values).
# -----------------------------
metabric_raw <- read.delim(metabric_expr_path, check.names = FALSE)
gene_col <- intersect(c("Hugo_Symbol", "GENE_SYMBOL", "gene_symbol"), colnames(metabric_raw))[1]
if (is.na(gene_col)) stop("Could not find a gene-symbol column in METABRIC expression file.")

known_id_cols <- c("Entrez_Gene_Id", "ENTREZ_GENE_ID", gene_col)
candidate_sample_cols <- setdiff(colnames(metabric_raw), known_id_cols)

suspicious_cols <- candidate_sample_cols[!grepl("^[A-Za-z0-9._-]+$", candidate_sample_cols)]
if (length(suspicious_cols) > 0) {
  stop(
    "METABRIC expression file contains column(s) that do not look like ",
    "plausible sample IDs and were NOT recognized as known ID columns: ",
    paste(suspicious_cols, collapse = ", "),
    ". Inspect the file before proceeding -- averaging an unexpected ",
    "metadata column into the expression matrix would silently corrupt ",
    "the analysis."
  )
}
cat("METABRIC expression file: ", length(candidate_sample_cols),
    " columns identified as sample expression columns after excluding known ID columns.\n", sep = "")

metabric_expr <- metabric_raw %>%
  dplyr::select(-any_of(known_id_cols[known_id_cols != gene_col])) %>%
  dplyr::rename(SYMBOL = all_of(gene_col)) %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

cat("METABRIC expression matrix:", dim(metabric_expr), "\n")

# -----------------------------
# Gene coverage check -- gated with explicit minimum
# -----------------------------
common_genes <- intersect(rownames(tcga_z), rownames(metabric_expr))
coverage_pct <- 100 * length(common_genes) / nrow(tcga_z)
cat("Gene coverage (TCGA top-5000 feature set found in METABRIC):",
    length(common_genes), "/", nrow(tcga_z), sprintf(" (%.1f%%)\n", coverage_pct))

MIN_COVERAGE_PCT <- 70
if (coverage_pct < MIN_COVERAGE_PCT) {
  stop(
    "Gene coverage (", round(coverage_pct, 1), "%) is below the pre-specified ",
    "minimum acceptable threshold (", MIN_COVERAGE_PCT, "%) for reliable ",
    "centroid-correlation projection."
  )
}

key_brca_genes <- c("ESR1", "PGR", "ERBB2", "MKI67", "GATA3", "FOXA1")
cat("\nKey BRCA gene presence in common gene set:\n")
print(setNames(key_brca_genes %in% common_genes, key_brca_genes))

# -----------------------------
# FIX: drop zero/NA-variance genes at the GENE level BEFORE scaling.
# scale() on a zero-variance row produces NaN for every sample in that
# row; filtering "complete" SAMPLES afterward (colSums(is.na(.))==0)
# would mark every sample incomplete if even one gene is degenerate.
# -----------------------------
metabric_common <- metabric_expr[common_genes, ]

gene_var <- apply(metabric_common, 1, var, na.rm = TRUE)
zero_var_genes <- names(gene_var)[is.na(gene_var) | gene_var == 0]
if (length(zero_var_genes) > 0) {
  cat(length(zero_var_genes), "gene(s) with zero/NA variance in METABRIC excluded",
      "before z-scoring (would otherwise silently drop every sample under a",
      "sample-wise NA filter).\n")
  metabric_common <- metabric_common[setdiff(rownames(metabric_common), zero_var_genes), ]
  common_genes <- setdiff(common_genes, zero_var_genes)
}

metabric_z <- t(scale(t(metabric_common)))

# Any remaining NAs at this point are gene-level (not sample-level)
# issues -- drop those specific genes, not samples.
bad_genes_post <- rownames(metabric_z)[rowSums(is.na(metabric_z)) > 0]
if (length(bad_genes_post) > 0) {
  cat(length(bad_genes_post), "additional gene(s) with residual NA after scaling excluded.\n")
  metabric_z <- metabric_z[setdiff(rownames(metabric_z), bad_genes_post), ]
  common_genes <- setdiff(common_genes, bad_genes_post)
}
cat("METABRIC z-scored matrix (post gene-level NA filter):", dim(metabric_z), "\n")
cat("Final common gene count after variance filtering:", length(common_genes), "\n")

# -----------------------------
# Nearest-centroid classification with gap-based confidence criterion
# -----------------------------
classify_to_centroids <- function(metabric_z, centroids, common_genes, label) {
  centroids_common <- centroids[common_genes, , drop = FALSE]
  cor_matrix <- cor(metabric_z[common_genes, ], centroids_common, method = "pearson")
  
  best_idx  <- max.col(cor_matrix, ties.method = "first")
  best_cor  <- cor_matrix[cbind(seq_len(nrow(cor_matrix)), best_idx)]
  second_best_cor <- apply(cor_matrix, 1, function(r) sort(r, decreasing = TRUE)[2])
  gap <- best_cor - second_best_cor
  
  assignment <- tibble(
    sample_id = rownames(cor_matrix),
    assigned_cluster = colnames(centroids_common)[best_idx],
    best_correlation = best_cor,
    second_best_correlation = second_best_cor,
    confidence_gap = gap
  )
  
  # Permutation null for the confidence gap: permute gene labels on the
  # METABRIC side, breaking true gene-to-gene correspondence with TCGA
  # centroids while preserving each sample's expression distribution.
  # NOTE: same offset (1300) is used for both k2 and k3 calls further
  # down, deliberately -- this means the SAME permuted gene orderings
  # are drawn both times, keeping the two null distributions built from
  # identical randomness rather than independent draws.
  set_pipeline_seed(offset = 1300)
  n_perm <- 200
  perm_gaps <- replicate(n_perm, {
    perm_genes <- sample(common_genes)
    cor_perm <- cor(metabric_z[perm_genes, ], centroids_common, method = "pearson")
    apply(cor_perm, 1, function(r) {
      sr <- sort(r, decreasing = TRUE); sr[1] - sr[2]
    })
  })
  null_gap_threshold <- quantile(as.vector(perm_gaps), 0.75)
  cat("[", label, "] Permutation-null confidence-gap threshold (75th pct): ",
      round(null_gap_threshold, 4), "\n", sep = "")
  
  assignment <- assignment %>%
    mutate(high_confidence = confidence_gap >= null_gap_threshold)
  
  cat("[", label, "] High-confidence assignments: ", sum(assignment$high_confidence),
      " / ", nrow(assignment), sprintf(" (%.1f%%)\n", 100 * mean(assignment$high_confidence)), sep = "")
  
  list(assignment = assignment, cor_matrix = cor_matrix, null_gap_threshold = null_gap_threshold)
}

class_k2 <- classify_to_centroids(metabric_z, centroids_k2, common_genes, "k2")
class_k3 <- classify_to_centroids(metabric_z, centroids_k3, common_genes, "k3")

write_csv(class_k2$assignment, here::here("results/tables/table_metabric_assignment_k2.csv"))
write_csv(class_k3$assignment, here::here("results/tables/table_metabric_assignment_k3.csv"))

p_gap_k3 <- ggplot(class_k3$assignment, aes(confidence_gap)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = class_k3$null_gap_threshold, linetype = "dashed", color = "red") +
  labs(title = "METABRIC classification confidence-gap distribution (k=3)",
       subtitle = "Dashed line: permutation-null-derived high-confidence threshold",
       x = "Confidence gap (top - second-best centroid correlation)") +
  theme_minimal()
ggsave(here::here("results/figures/figure_metabric_confidence_gap_k3.png"), p_gap_k3, width = 6.5, height = 5, dpi = 300)

# -----------------------------
# METABRIC clinical data + survival, matched covariate set
# -----------------------------
metabric_clinical <- read.delim(metabric_clin_path, comment.char = "#", check.names = FALSE)

# FIX: this duplicate-ID check now runs AFTER metabric_clinical exists
# (previously placed near the top of the script, before the object was
# created, which would throw "object not found" and stop the script).
if (any(duplicated(metabric_clinical$PATIENT_ID))) {
  warning(sum(duplicated(metabric_clinical$PATIENT_ID)),
          " duplicate PATIENT_ID(s) found in METABRIC clinical file -- ",
          "the join below may fan out. Investigate before trusting sample counts.")
}

unmatched_status <- metabric_clinical %>%
  filter(!is.na(OS_STATUS), !grepl("DECEASED|LIVING", OS_STATUS, ignore.case = TRUE)) %>%
  pull(OS_STATUS) %>% unique()
if (length(unmatched_status) > 0) {
  cat("WARNING: unmatched OS_STATUS values found (will become NA):",
      paste(unmatched_status, collapse = ", "), "\n")
}

metabric_clinical <- metabric_clinical %>%
  mutate(
    os_status_clean = case_when(
      grepl("DECEASED", OS_STATUS, ignore.case = TRUE) ~ 1,
      grepl("LIVING",   OS_STATUS, ignore.case = TRUE) ~ 0,
      TRUE ~ NA_real_
    ),
    os_months = suppressWarnings(as.numeric(OS_MONTHS)),
    os_years  = os_months / 12
  )
cat("\nos_status_clean breakdown (0=living, 1=deceased, NA=unresolved):\n")
print(table(metabric_clinical$os_status_clean, useNA = "always"))

# CRITICAL FIX: covariate columns are now chosen by DATA COMPLETENESS
# (pick_best_covariate(), R/utils.R), not by candidate-list order.
# FIX (this revision): er_status candidates now include "ER_IHC", the
# actual column name in this METABRIC clinical export. stage/grade
# candidates are left as-is for documentation purposes, but this export
# has NO stage or grade column at all (confirmed via colnames() check),
# so these will correctly resolve to NA and be excluded below -- this
# is a genuine data-availability gap, not a bug, and is documented in
# cross_platform_note above.
candidate_covariates <- c(
  age       = "AGE_AT_DIAGNOSIS",
  stage     = pick_best_covariate(c("TUMOR_STAGE", "STAGE"), metabric_clinical),
  grade     = pick_best_covariate(c("GRADE", "TUMOR_GRADE"), metabric_clinical),
  er_status = pick_best_covariate(c("ER_STATUS", "ER_IHC_STATUS", "ER_IHC"), metabric_clinical)
)
covariates_available <- candidate_covariates[!is.na(candidate_covariates)]
cat("\nMETABRIC covariates matched to TCGA discovery-cohort adjustment set\n",
    "(selected by data completeness, not list order):\n", sep = "")
print(covariates_available)
if (!("stage" %in% names(covariates_available)) || !("grade" %in% names(covariates_available))) {
  cat("\nNOTE: stage and/or grade unavailable in this METABRIC export -- the",
      "METABRIC multivariable Cox model below is NOT adjusted for the full",
      "age+stage set used in the TCGA discovery cohort. See",
      "NOTE_cross_platform_caveat.txt.\n")
}

# =============================================================
# run_metabric_survival_for_k()
#   cut_point: follow-up time (in os_years units) at which to split
#   follow-up for the time-split Cox sensitivity model. Defaults to 10
#   (approximate visual KM-curve crossing point for both k=2 and k=3 in
#   the current data) -- inspect figure_metabric_km_k*_all.png and
#   adjust per-k at the call site below if the crossing point differs.
# =============================================================
run_metabric_survival_for_k <- function(classification, k, cut_point = 10) {
  cat("\n========== METABRIC survival validation at k =", k, "==========\n")
  
  n_assignment <- nrow(classification$assignment)
  metabric_joined <- classification$assignment %>%
    inner_join(metabric_clinical, by = c("sample_id" = "PATIENT_ID"))
  n_after_join <- nrow(metabric_joined)
  
  metabric_cox_df <- metabric_joined %>%
    filter(!is.na(os_years), os_years >= 0, !is.na(os_status_clean)) %>%
    mutate(assigned_SNF_cluster = cluster_factor(assigned_cluster))
  n_after_filter <- nrow(metabric_cox_df)
  
  # Diagnostic block: verify filter-step row counts explicitly rather
  # than assuming a "zero exclusions" result is correct.
  cat("Row count diagnostics -- assignment:", n_assignment,
      "| after clinical join:", n_after_join,
      "| after survival-completeness filter:", n_after_filter, "\n")
  if (n_assignment != n_after_join) {
    warning("k=", k, ": row count changed during the clinical join (",
            n_assignment, " -> ", n_after_join, "). Check for duplicate ",
            "PATIENT_IDs (fan-out) or unmatched sample_ids (silent drops).")
  }
  cat("N with valid METABRIC survival data:", n_after_filter, "\n")
  
  fit_km <- survfit(Surv(os_years, os_status_clean) ~ assigned_SNF_cluster, data = metabric_cox_df)
  logrank <- survdiff(Surv(os_years, os_status_clean) ~ assigned_SNF_cluster, data = metabric_cox_df)
  
  logrank_p <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
  cat("Log-rank p (all assignments):", signif(logrank_p, 4), "\n")
  
  # FIX: KM plot now combined (curve + risk table) via ggpubr::ggarrange()
  # and saved with a file-size sanity check, rather than ggsave()'ing
  # p_km$plot alone (which either drops the risk table or, depending on
  # the graphics backend, can silently render a blank PNG -- observed in
  # script 10). k=3 also gets a wider canvas so its 3-cluster legend
  # isn't clipped.
  p_km <- ggsurvplot(fit_km, data = metabric_cox_df, pval = TRUE, risk.table = TRUE,
                     title = paste0("METABRIC OS by assigned SNF cluster (k=", k, ", all assignments)"))
  km_combined <- ggpubr::ggarrange(p_km$plot, p_km$table, ncol = 1, heights = c(3, 1))
  km_out_path <- here::here("results/figures", paste0("figure_metabric_km_k", k, "_all.png"))
  km_width <- if (k == 3) 9 else 7
  ggsave(km_out_path, plot = km_combined, width = km_width, height = 8, dpi = 300)
  stopifnot(file.exists(km_out_path), file.info(km_out_path)$size > 10000)
  cat("KM figure saved:", km_out_path, "(", file.info(km_out_path)$size, "bytes)\n")
  
  hc_df <- metabric_cox_df %>% filter(high_confidence)
  hc_df$assigned_SNF_cluster <- droplevels(hc_df$assigned_SNF_cluster)
  
  logrank_hc <- survdiff(Surv(os_years, os_status_clean) ~ assigned_SNF_cluster, data = hc_df)
  
  logrank_p_hc <- 1 - pchisq(logrank_hc$chisq, length(logrank_hc$n) - 1)
  cat("Log-rank p (high-confidence assignments only, N=", nrow(hc_df), "):",
      signif(logrank_p_hc, 4), "\n", sep = "")
  
  cox_formula_terms <- c("assigned_SNF_cluster")
  cox_df <- metabric_cox_df
  if ("age" %in% names(covariates_available)) {
    cox_df$age_covariate <- suppressWarnings(as.numeric(cox_df[[covariates_available["age"]]]))
    cox_formula_terms <- c(cox_formula_terms, "age_covariate")
  }
  clean_missing <- function(x) {
    x[x %in% c("", "NA", "[Not Available]", "Unknown", "N/A", "N/K")] <- NA
    x
  }
  
  if ("stage" %in% names(covariates_available)) {
    cox_df$stage_covariate <- factor(clean_missing(cox_df[[covariates_available["stage"]]]))
    cox_formula_terms <- c(cox_formula_terms, "stage_covariate")
  }
  if ("grade" %in% names(covariates_available)) {
    cox_df$grade_covariate <- factor(clean_missing(cox_df[[covariates_available["grade"]]]))
    cox_formula_terms <- c(cox_formula_terms, "grade_covariate")
  }
  if ("er_status" %in% names(covariates_available)) {
    cox_df$er_covariate <- factor(clean_missing(cox_df[[covariates_available["er_status"]]]))
    # FIX: relabel the "Positve" typo (raw cBioPortal field) at the
    # factor-level, not the raw vector -- this can't miss rows and
    # propagates automatically into broom::tidy() term names.
    levels(cox_df$er_covariate) <- dplyr::recode(levels(cox_df$er_covariate), "Positve" = "Positive")
    cox_formula_terms <- c(cox_formula_terms, "er_covariate")
    cat("ER covariate level counts:\n"); print(table(cox_df$er_covariate, useNA = "always"))
  }
  
  cox_df_complete <- cox_df %>%
    filter(if_all(any_of(c("age_covariate", "stage_covariate", "grade_covariate", "er_covariate")), ~ !is.na(.)))
  cat("Multivariable Cox N (complete cases on matched covariates):", nrow(cox_df_complete), "\n")
  
  cox_formula <- as.formula(paste0(
    "Surv(os_years, os_status_clean) ~ ", paste(cox_formula_terms, collapse = " + ")
  ))
  cox_multi <- coxph(cox_formula, data = cox_df_complete)
  cox_multi_tbl <- broom::tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(k = k, n = nrow(cox_df_complete), covariates_used = paste(names(covariates_available), collapse = "+"))
  
  ph_multi <- tryCatch(survival::cox.zph(cox_multi), error = function(e) NULL)
  if (!is.null(ph_multi)) {
    cat("\nProportional-hazards test (METABRIC multivariable Cox, k=", k, "):\n", sep = "")
    print(ph_multi)
  }
  ph_multi_tbl <- if (!is.null(ph_multi)) {
    as.data.frame(ph_multi$table) %>% rownames_to_column("term") %>% mutate(k = k)
  } else tibble()
  
  # ---------------------------------------------------------
  # NEW: Schoenfeld residual plot -- visualizes *how* the flagged
  # non-proportional effects (esp. assigned_SNF_cluster) drift over
  # follow-up time, rather than just reporting the cox.zph() p-value.
  # ---------------------------------------------------------
  if (!is.null(ph_multi)) {
    p_zph <- tryCatch(survminer::ggcoxzph(ph_multi), error = function(e) NULL)
    if (!is.null(p_zph)) {
      zph_out_path <- here::here("results/figures", paste0("figure_metabric_schoenfeld_k", k, ".png"))
      n_panels <- length(p_zph)
      png(zph_out_path, width = 7, height = 2.5 * n_panels, units = "in", res = 300)
      print(p_zph)
      dev.off()
      stopifnot(file.exists(zph_out_path), file.info(zph_out_path)$size > 10000)
      cat("Schoenfeld residual plot saved (", n_panels, "panel(s)):", zph_out_path, "\n")
    } else {
      cat("NOTE: ggcoxzph() failed for k=", k, " -- Schoenfeld plot not generated.\n", sep = "")
    }
  }
  
  # ---------------------------------------------------------
  # NEW: Time-split Cox model (cluster x period interaction) --
  # a PH-diagnostic sensitivity analysis. Follow-up is split at
  # cut_point (os_years); cluster is interacted with period so the HR
  # is allowed to differ before vs. after the split, rather than
  # forcing a single constant HR onto curves that are known to cross.
  # ---------------------------------------------------------
  cox_split_tbl <- tibble()
  split_df <- tryCatch(
    survival::survSplit(Surv(os_years, os_status_clean) ~ ., data = cox_df_complete,
                        cut = cut_point, episode = "period"),
    error = function(e) NULL
  )
  if (!is.null(split_df)) {
    split_df$period <- factor(split_df$period,
                              labels = c(paste0("<=", cut_point, "y"), paste0(">", cut_point, "y")))
    split_formula <- as.formula(paste0(
      "Surv(tstart, os_years, os_status_clean) ~ assigned_SNF_cluster * period + ",
      paste(setdiff(cox_formula_terms, "assigned_SNF_cluster"), collapse = " + ")
    ))
    cox_split <- tryCatch(coxph(split_formula, data = split_df), error = function(e) NULL)
    if (!is.null(cox_split)) {
      cox_split_tbl <- broom::tidy(cox_split, exponentiate = TRUE, conf.int = TRUE) %>%
        mutate(k = k, cut_point_years = cut_point, n = nrow(cox_df_complete))
      cat("\nTime-split Cox model (cluster x period interaction, cut at ", cut_point,
          "y), k=", k, ":\n", sep = "")
      print(cox_split_tbl)
    } else {
      cat("NOTE: time-split Cox model failed to fit for k=", k, ".\n", sep = "")
    }
  }
  
  # ---------------------------------------------------------
  # NEW: RMST pairwise comparisons -- a second, fully PH-assumption-free
  # sensitivity analysis. tau is set to the minimum of each cluster's
  # maximum observed follow-up time (the standard conservative default);
  # inspect for short-follow-up outliers pulling tau down if the result
  # looks truncated relative to the KM plots.
  # ---------------------------------------------------------
  tau <- floor(min(tapply(cox_df_complete$os_years, cox_df_complete$assigned_SNF_cluster, max)))
  cluster_levels <- levels(cox_df_complete$assigned_SNF_cluster)
  cluster_pairs <- combn(cluster_levels, 2, simplify = FALSE)
  
  rmst_results <- do.call(rbind, lapply(cluster_pairs, function(pair) {
    pair_df <- cox_df_complete[cox_df_complete$assigned_SNF_cluster %in% pair, ]
    arm <- as.numeric(pair_df$assigned_SNF_cluster == pair[2])
    out <- tryCatch(
      survRM2::rmst2(time = pair_df$os_years, status = pair_df$os_status_clean, arm = arm, tau = tau),
      error = function(e) NULL
    )
    if (is.null(out)) return(NULL)
    r <- out$unadjusted.result["RMST (arm=1)-(arm=0)", ]
    data.frame(k = k, group_ref = pair[1], group_comp = pair[2], tau_years = tau,
               rmst_diff_years = unname(r["Est."]),
               rmst_diff_lower = unname(r["lower .95"]),
               rmst_diff_upper = unname(r["upper .95"]),
               rmst_diff_p = unname(r["p"]))
  }))   # <-- lapply() AND do.call() both close here, cleanly
  
  # NOW rmst_results exists -- check it here, as its own statement
  if (any(is.na(rmst_results$rmst_diff_lower) | is.na(rmst_results$rmst_diff_upper))) {
    warning("k=", k, ": RMST confidence-interval bound(s) came back NA -- check column names in ",
            "survRM2::rmst2()'s output object against what this code indexes.")
  }
  
  cat("\nRMST pairwise comparisons (tau=", tau, " years), k=", k, ":\n", sep = "")
  print(rmst_results)
  
  list(
    logrank_p_all = logrank_p, n_all = nrow(metabric_cox_df),
    logrank_p_high_conf = logrank_p_hc, n_high_conf = nrow(hc_df),
    cox_multi = cox_multi_tbl, ph_multi = ph_multi_tbl,
    cox_split = cox_split_tbl, rmst = as_tibble(rmst_results)   # NEW
  )
}

surv_k2 <- run_metabric_survival_for_k(class_k2, 2, cut_point = 10)
surv_k3 <- run_metabric_survival_for_k(class_k3, 3, cut_point = 10)

external_validation_summary <- tibble(
  k = c(2, 3),
  n_total = c(surv_k2$n_all, surv_k3$n_all),
  logrank_p_all_assignments = c(surv_k2$logrank_p_all, surv_k3$logrank_p_all),
  n_high_confidence = c(surv_k2$n_high_conf, surv_k3$n_high_conf),
  logrank_p_high_confidence = c(surv_k2$logrank_p_high_conf, surv_k3$logrank_p_high_conf)
)
cat("\n\nExternal validation summary (k=2 vs k=3):\n"); print(external_validation_summary)
write_csv(external_validation_summary, here::here("results/tables/table_external_validation_summary_k2_k3.csv"))
write_csv(bind_rows(surv_k2$cox_multi, surv_k3$cox_multi),
          here::here("results/tables/table_metabric_cox_multivariable_k2_k3.csv"))
write_csv(bind_rows(surv_k2$ph_multi, surv_k3$ph_multi),
          here::here("results/tables/table_metabric_cox_ph_diagnostics_k2_k3.csv"))
write_csv(bind_rows(surv_k2$cox_split, surv_k3$cox_split),                          # NEW
          here::here("results/tables/table_metabric_cox_timesplit_k2_k3.csv"))      # NEW
write_csv(bind_rows(surv_k2$rmst, surv_k3$rmst),                                    # NEW
          here::here("results/tables/table_metabric_rmst_k2_k3.csv"))               # NEW

log_session_info(script_name, key_packages = c("survival", "survminer", "ggpubr", "survRM2"))
cat("\n\u2713 13_external_validation.R complete. Gene-symbol mapping fixed (0% -> real coverage);",
    "gene-level (not sample-level) variance filtering; KM figures fixed and size-checked;",
    "ER covariate matched (ER_IHC) and 'Positve' typo relabeled to 'Positive';",
    "Schoenfeld plots, time-split Cox, and RMST sensitivity analyses added to address",
    "confirmed non-proportional-hazards in the cluster term; stage/grade unavailability",
    "documented explicitly.\n")

close_logger(log_con, script_name)