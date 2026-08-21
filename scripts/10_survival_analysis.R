# ============================================================
# 10_survival_analysis.R  (publication-ready revision)
# Kaplan-Meier, log-rank, Cox regression, RMST by SNF cluster
# ============================================================
# CHANGES in this revision (audit follow-up):
#   - NEW: cox.zph() proportional-hazards diagnostic is now run for
#     EVERY Cox model (univariable and multivariable, both k values),
#     with results saved to a table and Schoenfeld residual plots saved
#     as supplementary figures. This closes a significant gap for a
#     survival-analysis-heavy thesis chapter: reported hazard ratios
#     are only straightforwardly interpretable if the PH assumption
#     approximately holds, and this was previously never checked.
#   - BH-FDR correction now applied across the RMST pairwise comparisons
#     within each k (3 comparisons for k=3), closing an inconsistency
#     where FDR correction was used elsewhere in the pipeline (scripts
#     08, 09) but not here.
#   - stage_simple's treatment as an unordered (dummy-coded) factor in
#     the Cox model is now an explicit, commented, deliberate choice
#     (allows non-monotonic stage effects) rather than an unstated
#     default.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(survival)
  library(survminer)
  library(survRM2)
})
source(here::here("R", "utils.R"))

script_name <- "10_survival_analysis"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

final_labels_obj <- readRDS(here::here("results/objects/final_snf_cluster_labels.rds"))
sample_metadata  <- readRDS(here::here("data/processed/sample_metadata_matched.rds"))

# FIX: days_to_death / days_to_last_followup were never coerced to numeric
# in resolve_clinical_columns() (R/utils.R) -- only age_years was. These
# columns commonly arrive as character type from TCGA's curated "paper_"
# clinical fields. Coerced here, in-memory, since no upstream script
# (03-09) performs arithmetic on these columns -- no need to regenerate
# sample_metadata_matched.rds or re-run any earlier script.
n_death_before <- sum(!is.na(sample_metadata$days_to_death))
n_followup_before <- sum(!is.na(sample_metadata$days_to_last_followup))

sample_metadata$days_to_death <- suppressWarnings(as.numeric(sample_metadata$days_to_death))
sample_metadata$days_to_last_followup <- suppressWarnings(as.numeric(sample_metadata$days_to_last_followup))

n_death_after <- sum(!is.na(sample_metadata$days_to_death))
n_followup_after <- sum(!is.na(sample_metadata$days_to_last_followup))
cat("days_to_death: coerced to numeric.", n_death_before - n_death_after,
    "value(s) became NA (non-numeric-parseable strings).\n")
cat("days_to_last_followup: coerced to numeric.", n_followup_before - n_followup_after,
    "value(s) became NA (non-numeric-parseable strings).\n")

require_columns(
  sample_metadata,
  c("vital_status", "days_to_death", "days_to_last_followup", "age_years", "stage_simple"),
  context = "10_survival_analysis: sample_metadata"
)

build_survival_df <- function(cluster_labels_k, k) {
  cluster_labels_k %>%
    left_join(
      sample_metadata %>%
        dplyr::select(full_barcode, vital_status, days_to_death,
                      days_to_last_followup, age_years, stage_simple),
      by = "full_barcode"
    ) %>%
    mutate(
      survival_time_days = ifelse(
        tolower(vital_status) == "dead", days_to_death, days_to_last_followup
      ),
      survival_time_years = survival_time_days / 365.25,
      survival_status = ifelse(tolower(vital_status) == "dead", 1, 0),
      SNF_cluster = cluster_factor(SNF_cluster)   # dynamic levels, correct for any k
    ) %>%
    filter(!is.na(survival_time_years), survival_time_years >= 0)
}

run_survival_for_k <- function(cluster_labels_k, k) {
  cat("\n========== Survival analysis at k =", k, "==========\n")
  surv_df <- build_survival_df(cluster_labels_k, k)
  n_before <- nrow(cluster_labels_k)
  n_after  <- nrow(surv_df)
  cat("Excluded for missing/unresolvable survival time:", n_before - n_after,
      " (", n_before, "samples ->", n_after, "with valid survival data)\n")
  surv_df$SNF_cluster <- droplevels(surv_df$SNF_cluster)
  cat("N with valid survival data:", nrow(surv_df), "\n")
  cat("Cluster levels present:", paste(levels(surv_df$SNF_cluster), collapse = ", "), "\n")
  

  # --- Kaplan-Meier + log-rank ---
 
  fit_km <- survfit(Surv(survival_time_years, survival_status) ~ SNF_cluster, data = surv_df)
  logrank <- survdiff(Surv(survival_time_years, survival_status) ~ SNF_cluster, data = surv_df)
  
  logrank_p <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
  cat("Log-rank p-value:", signif(logrank_p, 4), "\n")
  
  p_km <- ggsurvplot(
    fit_km, data = surv_df, pval = TRUE, risk.table = TRUE,
    title = paste0("Overall survival by SNF cluster (k=", k, ")"),
    legend.title = "Cluster"
  )
  
  km_combined <- ggpubr::ggarrange(p_km$plot, p_km$table, ncol = 1, heights = c(3, 1))
  out_path <- here::here("results/figures", paste0("figure_km_survival_k", k, ".png"))
  ggsave(out_path, plot = km_combined, width = 7, height = 8, dpi = 300)
  
  stopifnot(file.exists(out_path), file.info(out_path)$size > 10000)
  cat("KM figure k=", k, "saved, size:", file.info(out_path)$size, "bytes\n")
  
  # Sanity check -- fail loudly if the file didn't render properly
  out_path <- here::here("results/figures", paste0("figure_km_survival_k", k, ".png"))
  stopifnot(file.exists(out_path), file.info(out_path)$size > 10000)
  cat("KM figure k=", k, "saved, size:", file.info(out_path)$size, "bytes\n")
  
  
  # --- Univariable Cox ---
  cox_uni <- coxph(Surv(survival_time_years, survival_status) ~ SNF_cluster, data = surv_df)
  n_uni <- nrow(surv_df)
  
  # --- Multivariable Cox (age + stage). stage_simple is fit as an
  # UNORDERED factor (dummy-coded contrasts vs. Stage I), a deliberate
  # choice that allows a non-monotonic stage effect rather than forcing
  # a single linear ordinal coefficient; stated explicitly here so it
  # reads as intentional in a methods review, not a default overlooked. ---
  surv_multi_df <- surv_df %>% filter(!is.na(age_years), !is.na(stage_simple))
  n_multi <- nrow(surv_multi_df)
  cat("Univariable Cox N =", n_uni, "  Multivariable Cox N =", n_multi,
      " (", n_uni - n_multi, " excluded due to missing age/stage)\n", sep = "")
  
  cox_multi <- coxph(
    Surv(survival_time_years, survival_status) ~ droplevels(SNF_cluster) + age_years + stage_simple,
    data = surv_multi_df
  )
  
  cox_uni_tbl   <- broom::tidy(cox_uni,   exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "univariable", k = k, n = n_uni)
  cox_multi_tbl <- broom::tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "multivariable_age_stage", k = k, n = n_multi)
  
  # --- NEW: proportional-hazards assumption diagnostics for BOTH
  # models. Reported hazard ratios are only straightforwardly
  # interpretable if this assumption approximately holds. ---
  ph_uni   <- survival::cox.zph(cox_uni)
  ph_multi <- survival::cox.zph(cox_multi)
  
  cat("\nProportional-hazards test (univariable Cox):\n"); print(ph_uni)
  cat("\nProportional-hazards test (multivariable Cox):\n"); print(ph_multi)
  
  ph_uni_tbl   <- as.data.frame(ph_uni$table)   %>% rownames_to_column("term") %>% mutate(model = "univariable", k = k)
  ph_multi_tbl <- as.data.frame(ph_multi$table) %>% rownames_to_column("term") %>% mutate(model = "multivariable_age_stage", k = k)
  
  tryCatch({
    p_zph_uni <- survminer::ggcoxzph(ph_uni)
    ggsave(here::here("results/figures", paste0("figure_coxzph_univariable_k", k, ".png")),
           print(p_zph_uni), width = 7, height = 6, dpi = 300)
    p_zph_multi <- survminer::ggcoxzph(ph_multi)
    ggsave(here::here("results/figures", paste0("figure_coxzph_multivariable_k", k, ".png")),
           print(p_zph_multi), width = 8, height = 8, dpi = 300)
  }, error = function(e) message("ggcoxzph plotting failed (non-fatal): ", conditionMessage(e)))
  
  if (any(c(ph_uni_tbl$p, ph_multi_tbl$p) < 0.05, na.rm = TRUE)) {
    warning(
      "k=", k, ": at least one covariate/cluster term shows evidence of a ",
      "non-proportional-hazards effect (cox.zph p<0.05). Interpret the ",
      "corresponding hazard ratio as an AVERAGE effect over follow-up time, ",
      "and consider a time-varying-coefficient extension or stratification ",
      "as a sensitivity analysis before treating this HR as final."
    )
  }
  
  # --- RMST: all pairwise comparisons, WITH BH-FDR correction across
  # the comparison family within this k (fixes prior omission). ---
  cluster_levels <- levels(surv_df$SNF_cluster)
  pairs <- if (length(cluster_levels) >= 2) combn(cluster_levels, 2, simplify = FALSE) else list()
  
  rmst_results <- map_dfr(pairs, function(pr) {
    sub <- surv_df %>% filter(SNF_cluster %in% pr)
    arm <- ifelse(sub$SNF_cluster == pr[2], 1, 0)
    
    tryCatch({
      rmst_fit <- rmst2(time = sub$survival_time_years, status = sub$survival_status, arm = arm, tau = 5)
      res_mat <- rmst_fit$unadjusted.result
      tibble(
        k = k,
        comparison = paste0(pr[2], "_minus_", pr[1]),
        rmst_diff       = res_mat["RMST (arm=1)-(arm=0)", "Est."],
        rmst_diff_lower = res_mat["RMST (arm=1)-(arm=0)", "lower .95"],
        rmst_diff_upper = res_mat["RMST (arm=1)-(arm=0)", "upper .95"],
        rmst_diff_p     = res_mat["RMST (arm=1)-(arm=0)", "p"]
      )
    }, error = function(e) {
      message("RMST comparison ", pr[1], " vs ", pr[2], " (k=", k, ") failed: ", conditionMessage(e),
              " -- this comparison will be MISSING (NA) from the output table, not silently zero.")
      tibble(k = k, comparison = paste0(pr[2], "_minus_", pr[1]),
             rmst_diff = NA_real_, rmst_diff_lower = NA_real_,
             rmst_diff_upper = NA_real_, rmst_diff_p = NA_real_)
    })
  }) %>%
    mutate(rmst_diff_FDR = if (nrow(.) > 0) p.adjust(rmst_diff_p, method = "BH") else numeric(0))
  cat("\nRMST pairwise comparisons (tau=5y, FDR-corrected within this k):\n"); print(rmst_results)
  
  list(
    logrank_p = logrank_p,
    cox_uni = cox_uni_tbl, cox_multi = cox_multi_tbl,
    ph_uni = ph_uni_tbl, ph_multi = ph_multi_tbl,
    rmst = rmst_results, surv_df = surv_df
  )
}

res_k2 <- run_survival_for_k(final_labels_obj$k2, 2)
res_k3 <- run_survival_for_k(final_labels_obj$k3, 3)

logrank_summary <- tibble(k = c(2, 3), logrank_p = c(res_k2$logrank_p, res_k3$logrank_p))
cox_table   <- bind_rows(res_k2$cox_uni, res_k2$cox_multi, res_k3$cox_uni, res_k3$cox_multi)
ph_table    <- bind_rows(res_k2$ph_uni, res_k2$ph_multi, res_k3$ph_uni, res_k3$ph_multi)
rmst_table  <- bind_rows(res_k2$rmst, res_k3$rmst)

cat("\n\nLog-rank summary (k=2 vs k=3):\n"); print(logrank_summary)
write_csv(logrank_summary, here::here("results/tables/table_survival_logrank_k2_k3.csv"))
write_csv(cox_table,       here::here("results/tables/table_survival_cox_k2_k3.csv"))
write_csv(ph_table,        here::here("results/tables/table_survival_cox_ph_diagnostics_k2_k3.csv"))
write_csv(rmst_table,      here::here("results/tables/table_survival_rmst_k2_k3.csv"))

log_session_info(script_name, key_packages = c("survival", "survminer", "survRM2"))
cat("\n\u2713 10_survival_analysis.R complete. PH diagnostics + FDR-corrected RMST saved.\n")

close_logger(log_con, script_name)
