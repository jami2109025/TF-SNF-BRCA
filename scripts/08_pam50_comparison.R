# ============================================================
# 08_pam50_comparison.R  (publication-ready revision)
# Concordance between SNF clusters and PAM50 intrinsic subtypes
# ============================================================
# CHANGES in this revision (audit follow-up):
#   - Contingency table cell counts (and the expected-count minimum
#     used to decide chi-square vs. Fisher) are now explicitly saved,
#     so the test-selection logic is independently auditable rather
#     than only trusted from the branch behavior.
#   - Explicit interpretive guardrail persisted to disk: the permutation
#     empirical p-value demonstrates agreement BEYOND CHANCE between two
#     expression-derived classifiers, not validation against independent
#     ground truth (per the Script 02 caveat) -- this distinction is
#     easy to blur in a Results section and is now stated adjacent to
#     the numbers themselves, not only in a separate caveat file.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(mclust)
  library(aricode)
})
source(here::here("R", "utils.R"))

script_name <- "08_pam50_comparison"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

cat(readLines(here::here("results/tables/NOTE_pam50_independence_caveat.txt")), sep = "\n")
cat("\n")

final_labels_obj <- readRDS(here::here("results/objects/final_snf_cluster_labels.rds"))
sample_metadata  <- readRDS(here::here("data/processed/sample_metadata_matched.rds"))

run_pam50_comparison_for_k <- function(cluster_labels_k, k) {

  comparison_df <- cluster_labels_k %>%
    left_join(sample_metadata %>%
                dplyr::select(full_barcode, PAM50), by = "full_barcode") %>%
    filter(!is.na(PAM50), PAM50 != "")

  n_total <- nrow(cluster_labels_k)
  n_with_pam50 <- nrow(comparison_df)
  cat("\n--- k =", k, "---\n")
  cat("Samples with PAM50 available:", n_with_pam50, "of", n_total,
      sprintf(" (%.1f%%)", 100 * n_with_pam50 / n_total), "\n")

  contingency <- table(comparison_df$SNF_cluster, comparison_df$PAM50)
  cat("\nContingency table:\n"); print(contingency)

  chi_test <- suppressWarnings(chisq.test(contingency))
  min_expected <- min(chi_test$expected)
  use_fisher <- min_expected < 5
  cat("Minimum expected cell count:", round(min_expected, 2),
      " -> using", ifelse(use_fisher, "Fisher's exact (simulated)", "Chi-square"), "\n")

  set_pipeline_seed(offset = 800 + k)
  if (use_fisher) {
    test_result <- fisher.test(contingency, simulate.p.value = TRUE, B = 10000)
    test_used <- "Fisher's exact (Monte Carlo simulated p, B=10000, seeded)"
  } else {
    test_result <- chi_test
    test_used <- "Chi-square"
  }

  ari_value <- mclust::adjustedRandIndex(comparison_df$SNF_cluster, comparison_df$PAM50)
  nmi_value <- aricode::NMI(comparison_df$SNF_cluster, comparison_df$PAM50)

  set_pipeline_seed(offset = 900 + k)
  n_perm <- 1000
  perm_ari <- sapply(seq_len(n_perm), function(i) {
    mclust::adjustedRandIndex(comparison_df$SNF_cluster, sample(comparison_df$PAM50))
  })
  empirical_p <- mean(perm_ari >= ari_value)

  cat("ARI:", round(ari_value, 3), " NMI:", round(nmi_value, 3), "\n")
  cat("Permutation null (n=", n_perm, "): empirical p-value (ARI vs. chance) = ",
      ifelse(empirical_p == 0, paste0("<", 1 / n_perm), round(empirical_p, 4)), "\n", sep = "")
  cat("Association test used:", test_used, " p =", signif(test_result$p.value, 3), "\n")

  list(
    k = k, contingency = contingency, min_expected_cell = min_expected,
    test_used = test_used, p_value = test_result$p.value,
    ari = ari_value, nmi = nmi_value, empirical_p_ari = empirical_p,
    n_with_pam50 = n_with_pam50, n_total = n_total, comparison_df = comparison_df
  )
}

result_k2 <- run_pam50_comparison_for_k(final_labels_obj$k2, 2)
result_k3 <- run_pam50_comparison_for_k(final_labels_obj$k3, 3)

summary_table <- tibble(
  k = c(2, 3),
  n_with_pam50 = c(result_k2$n_with_pam50, result_k3$n_with_pam50),
  n_total = c(result_k2$n_total, result_k3$n_total),
  min_expected_cell_count = c(result_k2$min_expected_cell, result_k3$min_expected_cell),
  association_test = c(result_k2$test_used, result_k3$test_used),
  p_value = c(result_k2$p_value, result_k3$p_value),
  ARI = c(result_k2$ari, result_k3$ari),
  NMI = c(result_k2$nmi, result_k3$nmi),
  empirical_p_ARI_vs_permutation_null = c(result_k2$empirical_p_ari, result_k3$empirical_p_ari)
)
cat("\n\nSummary (k=2 vs k=3, both reported):\n"); print(summary_table)
write_csv(summary_table, here::here("results/tables/table_pam50_comparison_k2_k3.csv"))

write_csv(
  as.data.frame.matrix(result_k2$contingency) %>% rownames_to_column("SNF_cluster"),
  here::here("results/tables/table_pam50_contingency_k2.csv")
)
write_csv(
  as.data.frame.matrix(result_k3$contingency) %>% rownames_to_column("SNF_cluster"),
  here::here("results/tables/table_pam50_contingency_k3.csv")
)

writeLines(
  c("INTERPRETIVE GUARDRAIL for this script's results:",
    "===================================================",
    "The permutation empirical p-value demonstrates agreement BEYOND WHAT",
    "IS EXPECTED BY CHANCE between SNF clusters and PAM50 calls. Given the",
    "PAM50-non-independence caveat (see NOTE_pam50_independence_caveat.txt),",
    "this should be reported/discussed as: 'SNF clusters show significantly",
    "greater agreement with PAM50 subtype calls than expected by chance,'",
    "NOT as 'SNF clusters are validated against independent ground truth.'",
    "The two statements are not equivalent, and only the first is supported",
    "by this analysis."),
  here::here("results/tables/NOTE_pam50_result_interpretation.txt")
)

for (res in list(result_k2, result_k3)) {
  p <- ggplot(res$comparison_df, aes(x = SNF_cluster, fill = PAM50)) +
    geom_bar(position = "fill") +
    labs(title = paste0("PAM50 composition within SNF clusters (k=", res$k, ")"),
         y = "Proportion", x = "SNF cluster") +
    theme_minimal()
  ggsave(here::here("results/figures", paste0("figure_pam50_composition_k", res$k, ".png")),
         p, width = 6.5, height = 5, dpi = 300)
}

log_session_info(script_name, key_packages = c("mclust", "aricode"))
cat("\n\u2713 08_pam50_comparison.R complete.\n")

close_logger(log_con, script_name)
