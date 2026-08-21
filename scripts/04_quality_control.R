# ============================================================
# 04_quality_control.R  (publication-ready revision)
# QC: PCA, hierarchical clustering, outlier flagging, batch proxies
# ============================================================
# CHANGES in this revision (audit follow-up):
#   - Batch-proxy variables (possible_plate, possible_center) are now
#     written into sample_metadata_matched.rds itself (the object every
#     downstream script actually reads), NOT only into the previously
#     dead-end sample_metadata_with_qc.rds, so 06_snf_clustering.R can
#     formally test cluster assignment against them once cluster labels
#     exist (closing the audit's "batch proxies computed but never
#     tested against final clusters" gap).
#   - Tumor-purity overlay, when the optional file is supplied, is now
#     actually merged onto the PCA table and a purity-vs-PC1/PC2
#     correlation is reported, instead of only being loaded and saved
#     with no analysis performed on it.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})
source(here::here("R", "utils.R"))

script_name <- "04_quality_control"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

vst_top2000     <- readRDS(here::here("data/processed/vst_top2000_genes.rds"))
sample_metadata <- readRDS(here::here("data/processed/sample_metadata_matched.rds"))
stopifnot(identical(colnames(vst_top2000), rownames(sample_metadata)))

# -----------------------------
# PCA -- scaled and unscaled, both saved as an explicit sensitivity check
# -----------------------------
run_pca <- function(mat, scale_flag, label) {
  pca <- prcomp(t(mat), scale. = scale_flag)
  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
  pc_df <- as.data.frame(pca$x[, 1:5]) %>%
    rownames_to_column("full_barcode") %>%
    left_join(sample_metadata %>% dplyr::select(full_barcode, PAM50), by = "full_barcode")
  write_csv(pc_df, here::here("results/tables", paste0("table_pca_coords_", label, ".csv")))
  list(pca = pca, var_explained = var_explained, df = pc_df)
}

pca_scaled   <- run_pca(vst_top2000, scale_flag = TRUE,  label = "scaled")
pca_unscaled <- run_pca(vst_top2000, scale_flag = FALSE, label = "unscaled")

cat("Scaled PCA:   PC1 =", round(100 * pca_scaled$var_explained[1], 1),
    "% PC2 =", round(100 * pca_scaled$var_explained[2], 1), "%\n")
cat("Unscaled PCA: PC1 =", round(100 * pca_unscaled$var_explained[1], 1),
    "% PC2 =", round(100 * pca_unscaled$var_explained[2], 1), "%\n")

for (label_pca in list(list(d = pca_scaled, lab = "scaled"), list(d = pca_unscaled, lab = "unscaled"))) {
  p <- ggplot(label_pca$d$df, aes(PC1, PC2, color = PAM50)) +
    geom_point(alpha = 0.7, size = 1.8) +
    labs(
      title = paste0("PCA of top-2000-variance genes (", label_pca$lab, ")"),
      x = paste0("PC1 (", round(100 * label_pca$d$var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(100 * label_pca$d$var_explained[2], 1), "%)")
    ) +
    theme_minimal()
  ggsave(here::here("results/figures", paste0("figure_pca_", label_pca$lab, ".png")),
         p, width = 7, height = 5.5, dpi = 300)
}

# -----------------------------
# Hierarchical clustering (Ward.D2, Euclidean) with PAM50 annotation
# -----------------------------
dist_mat <- dist(t(vst_top2000), method = "euclidean")
hc       <- hclust(dist_mat, method = "ward.D2")

annotation_col <- sample_metadata %>%
  dplyr::select(PAM50) %>%
  mutate(PAM50 = ifelse(is.na(PAM50) | PAM50 == "", "Unknown", PAM50))

png(here::here("results/figures/figure_qc_hclust_samples.png"), width = 1400, height = 900, res = 150)
pheatmap(
  as.matrix(dist_mat), cluster_rows = hc, cluster_cols = hc,
  annotation_col = annotation_col, show_rownames = FALSE, show_colnames = FALSE,
  main = "Sample-sample Euclidean distance (Ward.D2), top-2000-variance genes"
)
dev.off()

# -----------------------------
# Outlier flagging (PCA z-score > 4) -- explicit retain/exclude decision
# -----------------------------
pc_z <- pca_scaled$df %>%
  mutate(z_PC1 = scale(PC1)[, 1], z_PC2 = scale(PC2)[, 1],
         possible_outlier = abs(z_PC1) > 4 | abs(z_PC2) > 4)

n_outliers <- sum(pc_z$possible_outlier)
cat("\nSamples flagged as possible PCA outliers (|z|>4):", n_outliers, "\n")
if (n_outliers > 0) print(pc_z %>% filter(possible_outlier) %>% dplyr::select(full_barcode, PAM50, z_PC1, z_PC2))

write_csv(pc_z, here::here("results/tables/table_qc_pca_outlier_flags.csv"))

decision_text <- c(
  paste0("QC outlier decision (", Sys.Date(), ")"),
  paste0("Samples flagged at |PCA z-score| > 4 on PC1/PC2: ", n_outliers),
  "DECISION: All samples retained in the primary analysis.",
  "RATIONALE: PCA outlier flags based on top-2000-variance gene expression",
  "can reflect genuine extreme molecular phenotypes rather than technical",
  "artifacts, and no independent technical QC metric was available to",
  "corroborate exclusion. Flagged sample IDs are reported in",
  "table_qc_pca_outlier_flags.csv. A sensitivity re-run of scripts 06-13",
  "excluding these samples IS PROVIDED as a supplementary analysis",
  "(see run_sensitivity_no_outliers.R) rather than only promised here."
)
writeLines(decision_text, here::here("results/tables/NOTE_qc_outlier_decision.txt"))

# -----------------------------
# Batch-proxy variables from TCGA barcode (correct field positions:
# Plate = chars 22-25, Center = chars 27-28)
# -----------------------------
sample_metadata$barcode_length <- nchar(sample_metadata$full_barcode)
cat("\nBarcode length distribution:\n"); print(table(sample_metadata$barcode_length))

full_length_ok <- sample_metadata$barcode_length >= 28

sample_metadata$possible_plate <- NA_character_
sample_metadata$possible_center <- NA_character_
sample_metadata$possible_plate[full_length_ok] <-
  substr(sample_metadata$full_barcode[full_length_ok], 22, 25)
sample_metadata$possible_center[full_length_ok] <-
  substr(sample_metadata$full_barcode[full_length_ok], 27, 28)

if (any(!full_length_ok)) {
  warning(sum(!full_length_ok), " sample barcode(s) shorter than the expected ",
          "28-character full TCGA aliquot barcode; plate/center left as NA for these.")
}

cat("\nPossible plate distribution (top 10):\n")
print(head(sort(table(sample_metadata$possible_plate), decreasing = TRUE), 10))
cat("\nPossible center distribution:\n")
print(table(sample_metadata$possible_center, useNA = "always"))

batch_df <- pca_scaled$df %>%
  left_join(sample_metadata %>%
              dplyr::select(full_barcode, possible_plate, possible_center),
            by = "full_barcode")

p_batch <- ggplot(batch_df, aes(PC1, PC2, color = possible_center)) +
  geom_point(alpha = 0.7, size = 1.6) +
  labs(title = "PCA colored by possible sequencing center (TCGA barcode chars 27-28)") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(here::here("results/figures/figure_qc_batch_by_center.png"), p_batch, width = 7, height = 6, dpi = 300)

# PC1/PC2 vs. batch-proxy association at the QC stage itself (an early
# warning even before clusters exist; the DEFINITIVE test against final
# cluster assignment is run in 06_snf_clustering.R once labels exist).
batch_pc_assoc <- test_cluster_batch_association(
  cluster_vector = cut(batch_df$PC1, breaks = 4),  # coarse PC1 bins as a quick proxy
  batch_vector   = batch_df$possible_center
)
cat("\nQuick early check -- PC1 (binned) vs. possible sequencing center:",
    batch_pc_assoc$test, " p =", signif(batch_pc_assoc$p_value, 3),
    "(a DEFINITIVE cluster-vs-batch test is run in 06_snf_clustering.R)\n")

# -----------------------------
# Tumor purity overlay -- now actually analyzed, not just loaded
# -----------------------------
purity_path <- here::here("data/external/tumor_purity_clean.csv")  # point at the CLEANED file
if (file.exists(purity_path)) {
  purity_df <- read_csv(purity_path, show_col_types = FALSE)
  cat("\nTumor purity file found and loaded:", nrow(purity_df), "rows.\n")
  
  purity_pc <- batch_df %>%
    mutate(sample_barcode_16 = substr(full_barcode, 1, 16)) %>%
    left_join(purity_df, by = "sample_barcode_16") %>%
    filter(!is.na(purity))
  
  cat("Samples matched to a purity value:", nrow(purity_pc), "of", nrow(batch_df), "\n")
  
  cor_pc1 <- suppressWarnings(cor.test(purity_pc$PC1, purity_pc$purity, method = "spearman"))
  cor_pc2 <- suppressWarnings(cor.test(purity_pc$PC2, purity_pc$purity, method = "spearman"))
  cat("Tumor purity vs PC1: rho =", round(cor_pc1$estimate, 3), " p =", signif(cor_pc1$p.value, 3), "\n")
  cat("Tumor purity vs PC2: rho =", round(cor_pc2$estimate, 3), " p =", signif(cor_pc2$p.value, 3), "\n")
  
  write_csv(purity_pc, here::here("results/tables/table_tumor_purity_merged.csv"))
} else {
  cat("\nNo tumor purity file found -- skipping purity QC.\n")
}

# Persist batch proxies onto the SAME metadata object every downstream
# script reads (sample_metadata_matched.rds), not only the previously
# dead-end sample_metadata_with_qc.rds -- this is the fix that makes
# the cluster-vs-batch test in 06_snf_clustering.R possible.
sample_metadata_matched <- readRDS(here::here("data/processed/sample_metadata_matched.rds"))
sample_metadata_matched$possible_plate  <- sample_metadata[rownames(sample_metadata_matched), "possible_plate"]
sample_metadata_matched$possible_center <- sample_metadata[rownames(sample_metadata_matched), "possible_center"]
saveRDS(sample_metadata_matched, here::here("data/processed/sample_metadata_matched.rds"))

saveRDS(sample_metadata, here::here("data/processed/sample_metadata_with_qc.rds"))
log_session_info(script_name)
cat("\n\u2713 04_quality_control.R complete. Batch proxies merged into sample_metadata_matched.rds.\n")

close_logger(log_con, script_name)
