# ============================================================
# 02_download_clinical_pam50.R  (publication-ready revision)
# Download TCGA clinical data and PAM50 subtype labels
# ============================================================
# CHANGES in this revision (audit follow-up):
#   - PAM50 duplicate-patient rows with CONFLICTING subtype calls are
#     now explicitly detected, counted, reported, and saved to a
#     dedicated table BEFORE deduplication, rather than silently
#     resolved by whichever row happens to appear first in the raw
#     TCGAquery_subtype() output order. Non-conflicting duplicates
#     (identical calls repeated) are still safely collapsed as before.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(TCGAbiolinks)
  library(tidyverse)
})
source(here::here("R", "utils.R"))

script_name <- "02_download_clinical_pam50"
log_con <- init_logger(script_name)
ensure_dirs()

# -----------------------------
# Clinical data
# -----------------------------
clinical_brca <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
cat("Clinical dimensions:", dim(clinical_brca), "\n")

clinical_clean <- data.frame(lapply(clinical_brca, function(x) {
  if (is.list(x)) sapply(x, function(v) paste(unlist(v), collapse = ";"))
  else x
}), check.names = FALSE)

write_csv(clinical_clean, here::here("data/raw/tcga_brca_clinical_raw.csv"))

# -----------------------------
# PAM50 subtype data
# -----------------------------
pam50_brca <- TCGAquery_subtype(tumor = "brca")
cat("\nPAM50 raw dimensions:", dim(pam50_brca), "\n")

pam50_patient_col_candidates <- grep(
  "barcode|patient|bcr", colnames(pam50_brca), ignore.case = TRUE, value = TRUE)
pam50_subtype_col_candidates <- grep(
  "subtype|pam50", colnames(pam50_brca), ignore.case = TRUE, value = TRUE)

if (length(pam50_patient_col_candidates) == 0 || length(pam50_subtype_col_candidates) == 0) {
  stop(
    "Could not auto-detect PAM50 patient-ID or subtype columns. ",
    "colnames(pam50_brca): ", paste(colnames(pam50_brca), collapse = ", ")
  )
}
pam50_patient_col <- pam50_patient_col_candidates[1]
pam50_subtype_col <- pam50_subtype_col_candidates[1]
cat("PAM50 patient column:", pam50_patient_col, "\n")
cat("PAM50 subtype column:", pam50_subtype_col, "\n")

pam50_brca$patient_id_12 <- substr(pam50_brca[[pam50_patient_col]], 1, 12)

# Explicitly detect and report CONFLICTING duplicate PAM50 calls
# (different subtype labels for the same patient) BEFORE collapsing,
# rather than letting distinct() silently keep an arbitrary row.
conflicts <- pam50_brca %>%
  group_by(patient_id_12) %>%
  filter(n_distinct(.data[[pam50_subtype_col]], na.rm = TRUE) > 1) %>%
  ungroup()

if (nrow(conflicts) > 0) {
  warning(
    n_distinct(conflicts$patient_id_12), " patient(s) have CONFLICTING PAM50 calls ",
    "across duplicate rows. These are saved to table_pam50_conflicting_duplicates.csv ",
    "for manual review; the FIRST row per patient is retained by default (documented, ",
    "not silent), but consider resolving these against the primary TCGA PAM50 ",
    "publication/source before finalizing."
  )
  write_csv(conflicts, here::here("results/tables/table_pam50_conflicting_duplicates.csv"))
} else {
  cat("No conflicting PAM50 duplicate calls detected.\n")
}

n_before <- nrow(pam50_brca)
pam50_dedup <- pam50_brca %>% distinct(patient_id_12, .keep_all = TRUE)
n_after <- nrow(pam50_dedup)
cat("PAM50 rows before dedup:", n_before, " after dedup:", n_after,
    " (", n_before - n_after, " duplicate patient rows removed, of which ",
    nrow(conflicts), " were label-conflicting)\n", sep = "")

write_csv(pam50_dedup, here::here("data/raw/tcga_brca_pam50_raw.csv"))

cat("\nPAM50 subtype distribution:\n")
print(table(pam50_dedup[[pam50_subtype_col]], useNA = "always"))

# -----------------------------
# Persistent methodological caveat (see audit, Script 02 + 08 + 13)
# -----------------------------
caveat_text <- c(
  "METHODOLOGICAL CAVEAT: PAM50 / CLAUDIN_SUBTYPE non-independence",
  "==================================================================",
  "PAM50 subtype calls (TCGAquery_subtype, used in 08_pam50_comparison.R)",
  "and METABRIC CLAUDIN_SUBTYPE calls (used in 13_external_validation.R)",
  "are both derived, by their respective classifiers, from tumor gene",
  "expression. They are NOT independent ground-truth labels relative to",
  "the SNF clusters produced by this pipeline, which are also derived",
  "from gene expression (plus TF activity, itself a function of",
  "expression). Concordance metrics (ARI, NMI, chi-square) between SNF",
  "clusters and PAM50/CLAUDIN_SUBTYPE should therefore be reported and",
  "discussed as AGREEMENT BETWEEN TWO EXPRESSION-BASED CLASSIFICATION",
  "SCHEMES, not as validation against an independent clinical or",
  "molecular ground truth. This caveat must be stated explicitly",
  "wherever such comparisons are reported in the manuscript."
)
writeLines(caveat_text, here::here("results/tables/NOTE_pam50_independence_caveat.txt"))

log_session_info(script_name, key_packages = "TCGAbiolinks")
cat("\n\u2713 02_download_clinical_pam50.R complete.\n")

close_logger(log_con, script_name)
