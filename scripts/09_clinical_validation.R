# ============================================================
# 09_clinical_validation.R  (publication-ready revision)
# Association between SNF clusters and clinicopathological variables
# ============================================================
# CHANGES in this revision (audit follow-up):
#   - Every test result row now reports N (valid, non-missing
#     observations), closing a silent gap where a heavily-missing
#     variable (e.g. `grade` in TCGA-BRCA, known to have substantial
#     missingness/legacy inconsistency) could produce a p-value with
#     no visible indication of how much data it was based on.
#   - test_categorical() now drops empty ROWS (zero-observation cluster
#     levels) in addition to the empty COLUMNS it already dropped, for
#     full symmetry/robustness (defensive; unlikely to trigger given
#     complete integer cluster labels, but costs nothing).
#   - FDR correction is explicitly documented as applied WITHIN each
#     k-solution's 4-test family (not pooled across k=2 and k=3), since
#     the two k-solutions are alternative partition hypotheses rather
#     than one combined testing family -- stated explicitly rather than
#     left as an implicit choice a reader would have to infer.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})
source(here::here("R", "utils.R"))

script_name <- "09_clinical_validation"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

final_labels_obj <- readRDS(here::here("results/objects/final_snf_cluster_labels.rds"))
sample_metadata  <- readRDS(here::here("data/processed/sample_metadata_matched.rds"))

require_columns(
  sample_metadata,
  c("stage_simple", "pathology_clean", "age_years", "grade"),
  context = "09_clinical_validation: sample_metadata"
)

test_categorical <- function(df, var, k, cluster_col = "SNF_cluster") {
  n_valid <- sum(!is.na(df[[var]]) & df[[var]] != "")
  tab <- table(df[[cluster_col]], df[[var]])
  tab <- tab[, colSums(tab) > 0, drop = FALSE]
  tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(tibble(variable = var, test = "insufficient data", statistic = NA, p_value = NA, n = n_valid))
  }
  chi <- suppressWarnings(chisq.test(tab))
  if (any(chi$expected < 5)) {
    set_pipeline_seed(offset = 950 + k)
    ft <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
    tibble(variable = var, test = "Fisher exact (simulated)", statistic = NA, p_value = ft$p.value, n = n_valid)
  } else {
    tibble(variable = var, test = "Chi-square", statistic = unname(chi$statistic), p_value = chi$p.value, n = n_valid)
  }
}

test_continuous <- function(df, var, cluster_col = "SNF_cluster") {
  d <- df %>% filter(!is.na(.data[[var]]))
  kw <- kruskal.test(d[[var]] ~ d[[cluster_col]])
  tibble(variable = var, test = "Kruskal-Wallis", statistic = unname(kw$statistic), p_value = kw$p.value, n = nrow(d))
}

run_clinical_validation_for_k <- function(cluster_labels_k, k) {
  cat("\n--- k =", k, "---\n")

  clinical_df <- cluster_labels_k %>%
    left_join(
      sample_metadata %>%
        dplyr::select(full_barcode, stage_simple, pathology_clean, age_years, grade),
      by = "full_barcode"
    )

  results <- bind_rows(
    test_categorical(clinical_df, "stage_simple", k),
    test_categorical(clinical_df, "pathology_clean", k),
    test_categorical(clinical_df, "grade", k),
    test_continuous(clinical_df, "age_years")
  ) %>%
    mutate(k = k, FDR = p.adjust(p_value, method = "BH"))
  print(results)
  list(results = results, clinical_df = clinical_df)
}

res_k2 <- run_clinical_validation_for_k(final_labels_obj$k2, 2)
res_k3 <- run_clinical_validation_for_k(final_labels_obj$k3, 3)

clinical_validation_table <- bind_rows(res_k2$results, res_k3$results)
cat("\n\nFull clinical validation table (k=2 and k=3):\n")
print(clinical_validation_table)
write_csv(clinical_validation_table, here::here("results/tables/table_clinical_validation_k2_k3.csv"))

writeLines(
  c("METHODS NOTE: multiple-testing correction scope.",
    "===================================================",
    "FDR (Benjamini-Hochberg) is applied WITHIN each k-solution's 4-test",
    "family (stage, pathology, grade, age), separately for k=2 and k=3.",
    "k=2 and k=3 are treated as ALTERNATIVE PARTITION HYPOTHESES, not a",
    "single pooled family of 8 comparisons, since they are co-equal",
    "parallel solutions rather than sequential/nested tests of the same",
    "hypothesis. This is a deliberate, stated methodological choice."),
  here::here("results/tables/NOTE_fdr_correction_scope.txt")
)

for (res in list(list(d = res_k2$clinical_df, k = 2), list(d = res_k3$clinical_df, k = 3))) {
  p_stage <- ggplot(res$d %>% filter(!is.na(stage_simple)), aes(x = SNF_cluster, fill = stage_simple)) +
    geom_bar(position = "fill") +
    labs(title = paste0("AJCC stage by SNF cluster (k=", res$k, ")"), y = "Proportion") +
    theme_minimal()
  ggsave(here::here("results/figures", paste0("figure_clinical_stage_k", res$k, ".png")),
         p_stage, width = 6.5, height = 5, dpi = 300)

  p_age <- ggplot(res$d %>% filter(!is.na(age_years)), aes(x = SNF_cluster, y = age_years)) +
    geom_boxplot(outlier.alpha = 0.4) +
    labs(title = paste0("Age at diagnosis by SNF cluster (k=", res$k, ")"), y = "Age (years)") +
    theme_minimal()
  ggsave(here::here("results/figures", paste0("figure_clinical_age_k", res$k, ".png")),
         p_age, width = 6, height = 5, dpi = 300)
}

log_session_info(script_name)
cat("\n\u2713 09_clinical_validation.R complete.\n")

close_logger(log_con, script_name)
