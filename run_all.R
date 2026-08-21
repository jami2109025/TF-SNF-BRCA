# ============================================================
# run_all.R
# Master driver: runs the full pipeline end-to-end in order, then
# performs the FULL-COVERAGE renv snapshot (see 00_setup.R's
# snapshot_full_pipeline() rationale -- a snapshot immediately after
# 00_setup.R alone would under-capture packages first used by later
# scripts).
# ============================================================

here::i_am("run_all.R")

scripts <- c(
  "scripts/00_setup.R",
  "scripts/01_download_tcga.R",
  "scripts/02_download_clinical_pam50.R",
  "scripts/03_preprocess_expression.R",
  "scripts/04_quality_control.R",
  "scripts/05_tf_activity.R",
  "scripts/06_snf_clustering.R",
  "scripts/07_benchmark_models.R",
  "scripts/08_pam50_comparison.R",
  "scripts/09_clinical_validation.R",
  "scripts/10_survival_analysis.R",
  "scripts/11_differential_expression_pathway.R",
  "scripts/12_master_regulator_network.R",
  "scripts/13_external_validation.R"
)

for (s in scripts) {
  cat("\n\n################################################################\n")
  cat("# Running:", s, "\n")
  cat("################################################################\n\n")
  t0 <- Sys.time()
  source(here::here(s), local = new.env(parent = globalenv()), echo = FALSE)
  cat("\n[", s, "] completed in",
      round(as.numeric(Sys.time() - t0, units = "mins"), 1), "minutes.\n")
}


cat("\n\nAll 13 scripts completed. Taking full-pipeline renv snapshot...\n")
source(here::here("scripts/00_setup.R"))  # brings snapshot_full_pipeline() into scope
snapshot_full_pipeline()

cat("\n\u2713 Full pipeline run complete. renv.lock now reflects every package used.\n")
