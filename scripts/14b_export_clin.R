# ============================================================
# 14b_export_clin.R
# Rebuild ONLY the METABRIC clinical frame from script 14 and export it
# for scripts/19_nested_bootstrap_sensitivity.R.
# ============================================================
# WHY THIS EXISTS
# ------------------------------------------------------------
# Script 19 needs script 14's prepared `clin` object. If your R session
# has been restarted since script 14 ran, `clin` is gone from memory.
#
# You do NOT need to re-run script 14 to get it back. Script 14's
# expensive stages are:
#
#   STEP 1  deriving TCGA cluster labels for 29 views  (~2 h, parLapply,
#           NO disk cache -- it would recompute in full)
#   STEP 3  loading and z-scoring the METABRIC expression matrix
#   STEP 5  the 1000-replicate paired bootstrap (~1 h; this one DOES
#           resume from results/objects/metabric_boot_*.rds)
#
# None of those touch `clin`. Its only inputs are the METABRIC clinical
# text file (~400 kB) and pick_best_covariate() from R/utils.R, so
# rebuilding it costs a few seconds.
#
# THE BLOCK BELOW IS COPIED VERBATIM FROM SCRIPT 14, STEP 4.
# Do not "improve" it. Its value is that it is byte-identical to the
# code that produced every METABRIC number in the manuscript, so the
# sensitivity analysis in script 19 cannot diverge from the primary
# analysis in script 14 on OS-status coding, on covariate selection, or
# on the cBioPortal "Positve" misspelling. If you ever edit STEP 4 in
# script 14, re-copy it here.
#
# HOW TO RUN
#   source(here::here("scripts/14b_export_clin.R"))
# Then run script 19.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})
source(here::here("R", "utils.R"))

metabric_clin_path <- here::here("data/external/metabric_clinical.txt")
if (!file.exists(metabric_clin_path)) {
  stop("METABRIC clinical file not found at: ", metabric_clin_path)
}

# ---- VERBATIM FROM 14_prognostic_added_value.R, STEP 4 --------------
metabric_clin <- read.delim(metabric_clin_path, check.names = FALSE, comment.char = "#")

id_col <- intersect(c("PATIENT_ID", "SAMPLE_ID", "Patient.ID"), colnames(metabric_clin))[1]
if (is.na(id_col)) stop("No patient/sample ID column found in METABRIC clinical file.")

os_months_col <- intersect(c("OS_MONTHS", "os_months"), colnames(metabric_clin))[1]
os_status_col <- intersect(c("OS_STATUS", "os_status"), colnames(metabric_clin))[1]
if (is.na(os_months_col) || is.na(os_status_col)) {
  stop("OS_MONTHS / OS_STATUS not found in METABRIC clinical file.")
}

age_col <- pick_best_covariate(c("AGE_AT_DIAGNOSIS", "Age.at.Diagnosis"), metabric_clin)
er_col  <- pick_best_covariate(c("ER_IHC", "ER_STATUS", "ER_status"), metabric_clin)
if (is.na(age_col) || is.na(er_col)) {
  stop("Could not locate an age and/or ER column in the METABRIC clinical file. ",
       "Adjusted analysis is the whole point of this script -- an unadjusted ",
       "comparison would reproduce the ER-driven result of script 13 and must ",
       "not be substituted silently. Columns available: ",
       paste(colnames(metabric_clin), collapse = ", "))
}
cat("METABRIC covariates selected by completeness -- age:", age_col, "| ER:", er_col, "\n")

clin <- metabric_clin %>%
  transmute(
    sample_id = .data[[id_col]],
    os_months = suppressWarnings(as.numeric(.data[[os_months_col]])),
    os_status_raw = as.character(.data[[os_status_col]]),
    age_covariate = suppressWarnings(as.numeric(.data[[age_col]])),
    er_covariate = as.character(.data[[er_col]])
  ) %>%
  mutate(
    os_status = case_when(
      grepl("^1|DECEASED|DEAD", os_status_raw, ignore.case = TRUE) ~ 1,
      grepl("^0|LIVING|ALIVE", os_status_raw, ignore.case = TRUE) ~ 0,
      TRUE ~ NA_real_
    ),
    os_years = os_months / 12,
    # cBioPortal's METABRIC export contains the misspelling "Positve";
    # left unhandled it silently creates a third ER level and changes
    # the model's degrees of freedom.
    er_covariate = case_when(
      grepl("^posit", er_covariate, ignore.case = TRUE) ~ "Positive",
      grepl("^negat", er_covariate, ignore.case = TRUE) ~ "Negative",
      TRUE ~ NA_character_
    )
  )

cat("METABRIC OS status resolved:\n"); print(table(clin$os_status, useNA = "ifany"))
# ---- END VERBATIM BLOCK --------------------------------------------

out_path <- here::here("results/objects/metabric_clin_prepared.rds")
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(clin, out_path)

# Fail loudly now rather than 8 hours into script 19's stage B.
required <- c("sample_id", "os_years", "os_status", "age_covariate", "er_covariate")
missing  <- setdiff(required, names(clin))
if (length(missing)) stop("Exported frame is missing columns: ",
                          paste(missing, collapse = ", "))

cat("\n=== EXPORT COMPLETE ===\n")
cat("Written:", out_path, "\n")
cat("Rows:", nrow(clin), "| columns:", paste(names(clin), collapse = ", "), "\n")
cat("Non-missing on all five model variables:",
    sum(stats::complete.cases(clin[, required])), "\n")
cat("ER levels:", paste(names(table(clin$er_covariate, useNA = "ifany")),
                        collapse = " / "), "\n")
cat("\nSanity check against script 14: after joining to the METABRIC\n")
cat("expression columns and filtering on complete outcome + covariates,\n")
cat("script 14 reported 1937 patients and 1125 events. Script 19 prints\n")
cat("the same two numbers at the start of its stage B -- if they differ,\n")
cat("stop and find out why before trusting any output.\n")
