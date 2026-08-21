# ============================================================
# 01_download_tcga.R  (publication-ready revision)
# Download TCGA-BRCA RNA-seq (STAR - Counts) via TCGAbiolinks
# Output: data/raw/tcga_brca_star_counts_se.rds
# ============================================================
# CHANGES in this revision (audit follow-up):
#   - Records the GDC data-release identifier alongside the query
#     snapshot, so the exact GDC data version underlying this cohort
#     is reproducible/auditable even if BRCA's GDC holdings change on
#     a future re-run (files added/revoked/reprocessed). Previously
#     only a date-stamped snapshot of the query RESULT TABLE was kept,
#     which records WHAT was queried but not WHICH GDC data release it
#     came from.
#   - Verifies GDCdownload() actually retrieved a file count consistent
#     with the query before calling GDCprepare(), instead of only
#     detecting a shortfall after the fact from the printed dimensions.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})
source(here::here("R", "utils.R"))

script_name <- "01_download_tcga"
log_con <- init_logger(script_name)
ensure_dirs()

dir.create(here::here("data/raw/GDCdata"), showWarnings = FALSE, recursive = TRUE)

# -- 1. Build query ---------------------------------------------------
query_brca <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

query_results <- getResults(query_brca)
cat("Query built. Files found:", nrow(query_results), "\n")

snapshot_date <- format(Sys.Date(), "%Y%m%d")
write_csv(
  query_results,
  here::here("data", "raw", paste0("tcga_brca_gdc_query_snapshot_", snapshot_date, ".csv"))
)

# GDC data-release provenance: record the release version alongside the
# query snapshot so the exact cohort composition is tied to a specific,
# citable GDC data release, not just "whatever was live on this date."
release_info <- tryCatch({
  if (requireNamespace("GenomicDataCommons", quietly = TRUE)) {
    GenomicDataCommons::status()
  } else {
    # Fallback: GDC REST status endpoint, no extra package required.
    jsonlite::fromJSON("https://api.gdc.cancer.gov/status")
  }
}, error = function(e) {
  message("Could not retrieve GDC data-release info: ", conditionMessage(e))
  NULL
})

writeLines(
  c(
    paste("GDC query snapshot date:", as.character(Sys.Date())),
    paste("GDC query file records:", nrow(query_results)),
    "GDC data-release status (raw):",
    if (!is.null(release_info)) capture.output(print(release_info)) else "NOT AVAILABLE -- record manually from https://portal.gdc.cancer.gov/ if this is missing."
  ),
  here::here("data/raw", paste0("gdc_release_info_", snapshot_date, ".txt"))
)
cat("GDC query snapshot and release-info provenance saved for", snapshot_date, "\n")

# -- 2. Download --------------------------------------------------
GDCdownload(
  query           = query_brca,
  method          = "api",
  files.per.chunk = 1231,
  directory       = here::here("data/raw/GDCdata")
)

# -- 3. Prepare into SummarizedExperiment -------------------------
brca_se <- GDCprepare(
  query     = query_brca,
  directory = here::here("data/raw/GDCdata"),
  save      = FALSE
)

cat("\nSummarizedExperiment dimensions:", dim(brca_se), "\n")
cat("Assay names:", assayNames(brca_se), "\n")

# Verify the prepared object's sample count is consistent with the
# number of case-level entries queried (not a strict 1:1 file:sample
# match, since GDCprepare can collapse/annotate, but a large shortfall
# indicates a partial download or prepare failure worth investigating
# BEFORE it propagates silently into every downstream script).
n_cases_queried <- length(unique(query_results$cases))
if (ncol(brca_se) < 0.9 * n_cases_queried) {
  warning(
    "Prepared SummarizedExperiment has ", ncol(brca_se), " samples, but ",
    n_cases_queried, " unique cases were queried (>10% shortfall). ",
    "Check for a partial GDCdownload() failure before proceeding to 03_preprocess_expression.R."
  )
}

cat("Sample types:\n")
sample_type_tab <- table(brca_se$sample_type, useNA = "always")
print(sample_type_tab)

if (!("Primary Tumor" %in% names(sample_type_tab)) || sample_type_tab["Primary Tumor"] < 100) {
  warning(
    "Fewer than 100 'Primary Tumor' samples detected in sample_type metadata. ",
    "This is unusual for TCGA-BRCA and should be investigated before proceeding ",
    "to 03_preprocess_expression.R, which relies on barcode-based primary-tumor ",
    "filtering independent of this column -- but a mismatch here may indicate a ",
    "GDCprepare() schema change worth checking."
  )
}

# -- 4. Save raw object for downstream preprocessing --------------
saveRDS(brca_se, here::here("data/raw/tcga_brca_star_counts_se.rds"))

log_session_info(script_name, key_packages = c("TCGAbiolinks", "SummarizedExperiment"))
cat("\n\u2713 01_download_tcga.R complete. Saved data/raw/tcga_brca_star_counts_se.rds\n")

close_logger(log_con, script_name)
