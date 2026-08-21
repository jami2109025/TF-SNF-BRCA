# ============================================================
# 17_crossover_confound_check.R   (NEW -- decisive test)
# Is the ten-year hazard crossover a CLUSTER finding, or a
# stromal/immune-content finding?
# ============================================================
# WHY THIS SCRIPT IS NOW THE MOST IMPORTANT ONE LEFT
# ------------------------------------------------------------
# Script 16 produced the thesis's one positive result: a hazard
# reversal at ten years that survives ER being given its OWN
# time-varying effect, corroborated model-free by RMST sign reversal
# (k=3: -1.08 years early, +3.24 years late).
#
# Script 15, run the same afternoon, produced the result that
# threatens it:
#
#     tumour purity by cluster, k=3:  Kruskal p = 2.6e-74
#         SNF_C1  median purity 0.818   Cohen's d = +1.24
#         SNF_C3  median purity 0.679   Cohen's d = -0.59
#     purity-residualised re-clustering, k=3:  ARI = 0.516
#
# A Cohen's d of 1.24 is a large effect, and an ARI of 0.516 means
# roughly half the k=3 partition dissolves once purity is regressed
# out. So the clusters are substantially a TUMOUR-CONTENT axis, and
# SNF_C3 -- precisely the cluster carrying the crossover -- is the
# LOW-purity, high-stromal/immune group.
#
# That raises an alternative explanation that a referee will reach for
# immediately, and which is entirely plausible biologically: immune
# infiltration is well known to associate with late outcome in breast
# cancer, so a low-purity cluster could show exactly this reversal
# without any novel subtype being involved.
#
# Two possible worlds, and they are reported very differently:
#
#   A. The crossover SURVIVES adjustment for stromal/immune content
#      -> the clusters carry time-varying prognostic information that
#         is not reducible to ER, not to existing subtype calls, and
#         not to tumour content. That is a genuine finding and the
#         centrepiece of the biology chapter.
#
#   B. The crossover does NOT survive
#      -> the finding is a stromal/immune-content effect. Still real
#         and still interesting, but it must be FRAMED that way, and
#         the thesis rests on the methodological contribution instead.
#
# Either answer is publishable. Claiming A without running this test is
# not, because the confound is already documented in this pipeline's
# own output.
#
# APPROACH
# METABRIC has no purity column, so tumour content is estimated from
# expression using Hallmark stromal and immune signatures scored per
# sample. The proxy is VALIDATED FIRST in TCGA, where true ABSOLUTE
# purity is available -- if the proxy does not correlate with measured
# purity there, it cannot be trusted in METABRIC and the script says so
# rather than proceeding.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(survival)
  library(org.Hs.eg.db)
})
source(here::here("R", "utils.R"))
source(here::here("R", "utils_benchmark.R"))

script_name <- "17_crossover_confound_check"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

CUT_YEARS <- 10
K_VALUES  <- c(2, 3)
MIN_PROXY_COR <- 0.30   # minimum |rho| in TCGA for the proxy to be usable

# ------------------------------------------------------------------
# STEP 1. Hallmark signatures for stromal and immune content
#
# Hallmark sets are used rather than the ESTIMATE gene lists because
# msigdbr is already a pipeline dependency and its version is captured
# by the archive in script 99, so the signature definition is frozen
# with everything else. Two stromal and two immune sets are averaged to
# avoid resting the adjustment on any single gene list.
# ------------------------------------------------------------------
get_hallmark <- function(set_names) {
  msig <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", category = "H"),
    error = function(e) msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  )
  sym_col <- intersect(c("gene_symbol", "human_gene_symbol"), colnames(msig))[1]
  set_col <- intersect(c("gs_name", "gs_id"), colnames(msig))[1]
  split(msig[[sym_col]], msig[[set_col]])[set_names]
}

STROMAL_SETS <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_ANGIOGENESIS")
IMMUNE_SETS  <- c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")

sigs <- get_hallmark(c(STROMAL_SETS, IMMUNE_SETS))
cat("Signature sizes:\n"); print(vapply(sigs, length, integer(1)))

# HARD CHECK. If msigdbr renames a set, or its schema changes, the
# lookup above returns NULL for that set. signature_score() would then
# quietly return all-NA, score_block() would propagate NA, every patient
# would be dropped by the complete-case filter, and the script could
# still print a VERDICT computed from whatever survived. A silent
# adjustment-for-nothing is exactly the failure mode that would
# manufacture a false "the crossover survives". Fail loudly instead.
missing_sets <- names(sigs)[vapply(sigs, function(x) is.null(x) || length(x) < 20, logical(1))]
if (length(missing_sets) > 0 || length(sigs) != length(c(STROMAL_SETS, IMMUNE_SETS))) {
  stop("Hallmark signature(s) missing or too small: ",
       paste(c(missing_sets, "(check msigdbr set names)"), collapse = ", "),
       ". Available Hallmark sets can be listed with: ",
       "unique(msigdbr::msigdbr(species='Homo sapiens', category='H')$gs_name)")
}

#' Mean per-gene z-score across a signature, computed WITHIN cohort.
#' Deliberately simple: a mean of z-scores is transparent, has no tuning
#' parameters, and cannot silently fail the way a rank-based enrichment
#' method can on a platform it was not designed for. The claim being
#' tested is a coarse one (is there residual stromal/immune signal?), so
#' a coarse, robust score is appropriate.
signature_score <- function(expr_mat, genes) {
  g <- intersect(genes, rownames(expr_mat))
  if (length(g) < 20) return(rep(NA_real_, ncol(expr_mat)))
  z <- t(scale(t(expr_mat[g, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

score_block <- function(expr_mat, label = "") {
  strom <- rowMeans(sapply(STROMAL_SETS, function(s) signature_score(expr_mat, sigs[[s]])), na.rm = TRUE)
  immun <- rowMeans(sapply(IMMUNE_SETS,  function(s) signature_score(expr_mat, sigs[[s]])), na.rm = TRUE)
  frac_na <- mean(!is.finite(strom) | !is.finite(immun))
  if (frac_na > 0.05) {
    stop("score_block(", label, "): ", round(100 * frac_na, 1),
         "% of samples have a non-finite stromal or immune score. The",
         " signature genes are probably poorly represented in this",
         " expression matrix. Do not proceed -- adjusting for a degenerate",
         " score is equivalent to not adjusting at all.")
  }
  tibble(sample_id = colnames(expr_mat), stromal_score = strom, immune_score = immun)
}

# ------------------------------------------------------------------
# STEP 2. VALIDATE THE PROXY IN TCGA, where real purity is available.
#
# This gate exists because an unvalidated proxy that fails to capture
# tumour content would produce a FALSE "the crossover survives
# adjustment" result -- adjusting for noise adjusts for nothing. The
# script refuses to proceed if the proxy does not track measured purity.
# ------------------------------------------------------------------
cat("\n=== STEP 2: validating the stromal/immune proxy against ABSOLUTE purity in TCGA ===\n")

vst_all <- readRDS(here::here("data/processed/vst_matrix_all_genes.rds"))
mapped  <- map_to_symbol_dedup(vst_all, variance_reference = vst_all, org_db = org.Hs.eg.db)
tcga_sym <- mapped$matrix

purity_path <- here::here("results/tables/table_tumor_purity_merged.csv")
if (!file.exists(purity_path)) stop("table_tumor_purity_merged.csv not found; run script 04/15 first.")
purity_df <- read_csv(purity_path, show_col_types = FALSE) %>%
  dplyr::select(full_barcode, purity) %>% filter(!is.na(purity))

tcga_scores <- score_block(tcga_sym, "TCGA") %>%
  inner_join(purity_df, by = c("sample_id" = "full_barcode"))

rho_strom <- cor(tcga_scores$stromal_score, tcga_scores$purity, method = "spearman", use = "complete.obs")
rho_immun <- cor(tcga_scores$immune_score,  tcga_scores$purity, method = "spearman", use = "complete.obs")
cat(sprintf("  n = %d\n  stromal score vs purity: rho = %+.3f\n  immune  score vs purity: rho = %+.3f\n",
            nrow(tcga_scores), rho_strom, rho_immun))
cat("  (negative rho is expected: high stromal/immune content means LOW tumour purity)\n")

proxy_ok <- (abs(rho_strom) >= MIN_PROXY_COR) || (abs(rho_immun) >= MIN_PROXY_COR)
if (!proxy_ok) {
  stop("The stromal/immune proxy does not track measured purity in TCGA ",
       "(|rho| < ", MIN_PROXY_COR, "). Adjusting METABRIC for this score would ",
       "be adjusting for noise and could produce a spurious 'the crossover ",
       "survives' result. Do not proceed; use a validated deconvolution ",
       "method instead, or report the confound as an unresolved limitation.")
}
cat("  -> proxy VALIDATED; safe to use in METABRIC.\n")

write_csv(tibble(score = c("stromal", "immune"), spearman_vs_purity = c(rho_strom, rho_immun),
                 n = nrow(tcga_scores)),
          here::here("results/tables/table_stromal_immune_proxy_validation.csv"))

# ------------------------------------------------------------------
# STEP 3. Score METABRIC
# ------------------------------------------------------------------
cat("\n=== STEP 3: scoring METABRIC ===\n")
metabric_raw <- read.delim(here::here("data/external/metabric_expression.txt"), check.names = FALSE)
gene_col <- intersect(c("Hugo_Symbol", "GENE_SYMBOL", "gene_symbol"), colnames(metabric_raw))[1]
metabric_expr <- metabric_raw %>%
  dplyr::select(-any_of(setdiff(c("Entrez_Gene_Id", "ENTREZ_GENE_ID", gene_col), gene_col))) %>%
  dplyr::rename(SYMBOL = all_of(gene_col)) %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("SYMBOL") %>% as.matrix()

mb_scores <- score_block(metabric_expr, "METABRIC")
cat("  METABRIC scored:", nrow(mb_scores), "samples\n")

# ------------------------------------------------------------------
# STEP 4. Clinical + assignments (same construction as script 16)
# ------------------------------------------------------------------
metabric_clin <- read.delim(here::here("data/external/metabric_clinical.txt"),
                            check.names = FALSE, comment.char = "#")
id_col <- intersect(c("PATIENT_ID", "SAMPLE_ID"), colnames(metabric_clin))[1]

clin <- metabric_clin %>%
  transmute(
    sample_id = .data[[id_col]],
    os_years  = suppressWarnings(as.numeric(.data[["OS_MONTHS"]])) / 12,
    os_raw    = as.character(.data[["OS_STATUS"]]),
    age       = suppressWarnings(as.numeric(.data[["AGE_AT_DIAGNOSIS"]])),
    er_raw    = as.character(.data[["ER_IHC"]]),
    subtype   = as.character(.data[["CLAUDIN_SUBTYPE"]])
  ) %>%
  mutate(
    os_status = case_when(grepl("^1|DECEASED|DEAD", os_raw, ignore.case = TRUE) ~ 1,
                          grepl("^0|LIVING|ALIVE", os_raw, ignore.case = TRUE) ~ 0,
                          TRUE ~ NA_real_),
    er = case_when(grepl("^posit", er_raw, ignore.case = TRUE) ~ "Positive",
                   grepl("^negat", er_raw, ignore.case = TRUE) ~ "Negative",
                   TRUE ~ NA_character_),
    subtype = ifelse(is.na(subtype) | subtype %in% c("", "NC"), NA_character_, subtype)
  )

results <- list()

for (k in K_VALUES) {
  asg <- read_csv(here::here(sprintf("results/tables/table_metabric_assignment_k%d.csv", k)),
                  show_col_types = FALSE) %>%
    dplyr::select(sample_id, assigned_cluster)

  d <- clin %>%
    inner_join(asg, by = "sample_id") %>%
    inner_join(mb_scores, by = "sample_id") %>%
    filter(is.finite(os_years), os_years > 0, !is.na(os_status), !is.na(age), !is.na(er),
           !is.na(subtype), !is.na(stromal_score), !is.na(immune_score)) %>%
    mutate(cl = factor(assigned_cluster), er = factor(er), subtype = factor(subtype),
           stromal_z = as.numeric(scale(stromal_score)),
           immune_z  = as.numeric(scale(immune_score))) %>%
    dplyr::select(-os_raw, -er_raw)

  cat("\n========== k =", k, "==========\n")
  cat("N:", nrow(d), " events:", sum(d$os_status), "\n")

  # Do the clusters differ in the proxy, as they do in real purity in TCGA?
  kw_s <- kruskal.test(stromal_score ~ cl, data = d)
  kw_i <- kruskal.test(immune_score  ~ cl, data = d)
  cat(sprintf("Cluster differences in proxy -- stromal p = %.3g | immune p = %.3g\n",
              kw_s$p.value, kw_i$p.value))
  print(d %>% group_by(cl) %>%
          summarise(n = n(), median_stromal = median(stromal_score),
                    median_immune = median(immune_score), .groups = "drop"))

  dsplit <- survSplit(Surv(os_years, os_status) ~ ., data = d,
                      cut = CUT_YEARS, episode = "period_id", id = "subject_id") %>%
    mutate(period = factor(ifelse(period_id == 1, "early_0_10y", "late_gt10y"),
                           levels = c("early_0_10y", "late_gt10y")))

  # T4  = script 16's model: cluster x period + ER x period + age
  # T5  = T4 + stromal x period + immune x period
  #
  # The proxy terms are given their OWN period interactions, not just
  # main effects. Adjusting only for a main effect would leave a
  # time-varying stromal effect free to load onto the cluster:period
  # term -- the identical mistake that fitting T1 instead of T4 would
  # have made with ER.
  T4 <- coxph(Surv(tstart, os_years, os_status) ~ cl * period + er * period + age, data = dsplit)
  T5 <- coxph(Surv(tstart, os_years, os_status) ~ cl * period + er * period + age +
                stromal_z * period + immune_z * period, data = dsplit)

  # T6, the strictest model available: every known competing explanation
  # is given its own time-varying effect simultaneously -- ER, claudin
  # subtype, stromal content and immune content. Script 16 established
  # that subtype carries real prognostic information (M2 vs M1,
  # p = 8.8e-6), so leaving it out of the adjusted model would leave an
  # obvious gap. A cluster x period interaction that survives T6 is not
  # explainable by anything this dataset can measure.
  T6 <- coxph(Surv(tstart, os_years, os_status) ~ cl * period + er * period + age +
                stromal_z * period + immune_z * period + subtype * period, data = dsplit)

  tidy_keep <- function(fit, nm) {
    broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(model = nm, k = k, .before = 1)
  }
  both <- bind_rows(tidy_keep(T4, "T4_no_content_adjustment"),
                    tidy_keep(T5, "T5_adjusted_for_stromal_immune"),
                    tidy_keep(T6, "T6_adjusted_for_content_and_subtype"))
  results[[as.character(k)]] <- both

  cat("\nT5 (cluster x period, with stromal and immune ALSO time-varying):\n")
  print(T5 %>% broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
          dplyr::select(term, estimate, conf.low, conf.high, p.value), n = 30)

  keyterm <- function(mod) both %>% filter(model == mod, grepl("^cl.*:period", term))
  key4 <- keyterm("T4_no_content_adjustment")
  key5 <- keyterm("T5_adjusted_for_stromal_immune")
  key6 <- keyterm("T6_adjusted_for_content_and_subtype")

  cat("\nCluster x period interaction across the adjustment ladder:\n")
  print(bind_rows(
    key4 %>% mutate(adjustment = "T4: ER"),
    key5 %>% mutate(adjustment = "T5: ER + stromal + immune"),
    key6 %>% mutate(adjustment = "T6: ER + stromal + immune + subtype")
  ) %>% dplyr::select(adjustment, term, estimate, conf.low, conf.high, p.value), n = 30)

  # An EMPTY key table means the interaction term was not estimated at
  # all -- a model failure, not a negative result. Reporting "does not
  # survive" in that case would be a false negative, so it is separated
  # out explicitly.
  if (nrow(key5) == 0 || nrow(key6) == 0) {
    cat("\n>>> k =", k, ": FAILED -- no cluster x period coefficient was estimated.",
        "Inspect the model; do NOT read this as a negative result.\n")
  } else {
    surv5 <- any(key5$p.value < 0.05, na.rm = TRUE)
    surv6 <- any(key6$p.value < 0.05, na.rm = TRUE)
    which6 <- key6$term[which(key6$p.value < 0.05)]

    cat("\n>>> VERDICT k = ", k, "\n", sep = "")
    cat("    T5 (ER + stromal + immune)            : ",
        if (surv5) "SURVIVES" else "does NOT survive", "\n", sep = "")
    cat("    T6 (+ claudin subtype, strictest)     : ",
        if (surv6) "SURVIVES" else "does NOT survive", "\n", sep = "")
    if (surv6) {
      cat("    surviving term(s): ", paste(which6, collapse = ", "), "\n", sep = "")
      cat("    -> Not reducible to ER, to existing subtype calls, or to tumour\n")
      cat("       content. Report it as the biological finding.\n")
    } else if (surv5) {
      cat("    -> Survives tumour content but NOT the addition of subtype. The\n")
      cat("       crossover is shared with existing subtype calls; report it as\n")
      cat("       consistent with them, not as independent of them.\n")
    } else {
      cat("    -> Explained by stromal/immune content. Frame it as a\n")
      cat("       tumour-microenvironment effect, NOT a novel subtype property.\n")
    }
  }
}

all_results <- bind_rows(results)
write_csv(all_results, here::here("results/tables/table_crossover_content_adjustment.csv"))

writeLines(c(
  "IS THE TEN-YEAR CROSSOVER A CLUSTER FINDING OR A TUMOUR-CONTENT FINDING?",
  "======================================================================",
  "",
  "WHY THIS TEST WAS NECESSARY.",
  "Script 16 found a hazard reversal at ten years that survives ER being",
  "given its own time-varying effect, corroborated model-free by RMST",
  "sign reversal. Script 15 found that the same clusters differ markedly",
  "in tumour purity (k=3 Kruskal p = 2.6e-74; SNF_C1 Cohen's d = +1.24)",
  "and that the k=3 partition only partly survives purity residualisation",
  "(ARI = 0.516). The cluster carrying the crossover, SNF_C3, is the",
  "LOW-purity, high-stromal/immune group. Immune infiltration is already",
  "known to associate with late outcome in breast cancer, so a referee",
  "will propose that explanation immediately -- and would be right to.",
  "",
  "WHAT WAS DONE.",
  "METABRIC has no purity column, so stromal and immune content were",
  "scored from expression using Hallmark signatures (EMT + angiogenesis;",
  "inflammatory response + interferon-gamma response) as the mean of",
  "per-gene z-scores. The proxy was VALIDATED FIRST in TCGA against",
  "measured ABSOLUTE purity -- see",
  "table_stromal_immune_proxy_validation.csv. The script refuses to",
  "proceed if the proxy fails that check, because adjusting for a score",
  "that does not capture tumour content would adjust for noise and could",
  "manufacture a false 'the crossover survives' result.",
  "",
  "The adjusted model (T5) gives the stromal and immune scores their OWN",
  "period interactions, not merely main effects. Adjusting only for main",
  "effects would leave a time-varying stromal effect free to load onto",
  "the cluster x period term -- exactly the error that fitting T1 rather",
  "than T4 would have made with ER.",
  "",
  "HOW TO REPORT EITHER OUTCOME.",
  "SURVIVES: the clusters carry time-varying prognostic information not",
  "reducible to ER, to existing subtype calls, or to tumour content.",
  "This is the biology chapter, and the purity association becomes a",
  "characterisation of the clusters rather than a threat to them.",
  "",
  "DOES NOT SURVIVE: the crossover is a tumour-microenvironment effect.",
  "Report it as such -- it remains a real and interesting observation,",
  "but it is not a novel subtype property, and the thesis rests on the",
  "methodological contribution (the equivalence-tested added-value",
  "framework) instead.",
  "",
  "LIMITATION. A signature-based score is a coarse proxy for cell-type",
  "composition. A full deconvolution (CIBERSORTx, ESTIMATE, xCell) would",
  "be more precise and is the obvious extension. The TCGA validation",
  "step bounds how badly the proxy can be failing, but does not make it",
  "exact."
), here::here("results/tables/NOTE_crossover_confound_check.txt"))

log_session_info(script_name, key_packages = c("survival", "msigdbr", "broom"))
cat("\n✓ 17_crossover_confound_check.R complete.\n")

close_logger(log_con, script_name)
