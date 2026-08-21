# ============================================================
# 03_preprocess_expression.R  (publication-ready revision)
# Preprocess TCGA-BRCA expression data
# ============================================================
# CHANGES in this revision (audit follow-up):
#   - Counts are now rounded to integers BEFORE the low-count gene
#     filter is applied (previously the filter ran on unrounded STAR
#     multi-mapped fractional counts, and rounding happened only
#     afterward for DESeq2 -- a second-order but avoidable
#     inconsistency between the filter and the model that consumes it).
#   - Low-count filter minimum-sample threshold now uses ceiling()
#     instead of round(), so the "≥10% of samples" criterion cannot be
#     silently relaxed by round-half-down behavior.
#   - Explicit unsupervised-discovery framing note persisted to disk:
#     feature selection (top2000/top5000 variance genes), clustering,
#     and downstream marker-gene testing all draw on the SAME cohort
#     with no internal holdout. This is standard and defensible for a
#     discovery study, but must be stated as a design property, not
#     left implicit, since it is the root of the circularity caveats
#     that scripts 06/08/11/12 each document downstream.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(SummarizedExperiment)
  library(DESeq2)
  library(tidyverse)
})
source(here::here("R", "utils.R"))

script_name <- "03_preprocess_expression"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

# -----------------------------
# Load data
# -----------------------------
brca_se       <- readRDS(here::here("data/raw/tcga_brca_star_counts_se.rds"))
clinical_brca <- read.csv(here::here("data/raw/tcga_brca_clinical_raw.csv"))
pam50_brca    <- read.csv(here::here("data/raw/tcga_brca_pam50_raw.csv"))

# -----------------------------
# Extract count matrix
# -----------------------------
if ("unstranded" %in% assayNames(brca_se)) {
  count_matrix <- assay(brca_se, "unstranded")
} else {
  count_matrix <- assay(brca_se, 1)
}
count_matrix    <- as.matrix(count_matrix)
sample_metadata <- as.data.frame(colData(brca_se))
gene_annotation <- as.data.frame(rowData(brca_se))

preprocess_summary <- tibble(
  step = "Raw download", genes = nrow(count_matrix), samples = ncol(count_matrix)
)

# -----------------------------
# Primary tumor filter (barcode pos 14-15 == "01")
# -----------------------------
sample_type_code <- substr(colnames(count_matrix), 14, 15)
cat("Sample type codes:\n"); print(table(sample_type_code))

keep_primary    <- sample_type_code == "01"
count_matrix    <- count_matrix[, keep_primary]
sample_metadata <- sample_metadata[keep_primary, ]

preprocess_summary <- bind_rows(preprocess_summary, tibble(
  step = "Primary tumor only", genes = nrow(count_matrix), samples = ncol(count_matrix)
))

# -----------------------------
# Remove duplicate-patient aliquots (keep highest library size)
# -----------------------------
duplicate_table <- tibble(
  full_barcode = colnames(count_matrix),
  patient_id   = substr(colnames(count_matrix), 1, 12),
  library_size = colSums(count_matrix)
)

samples_to_keep <- duplicate_table %>%
  group_by(patient_id) %>%
  slice_max(order_by = library_size, n = 1, with_ties = FALSE) %>%
  pull(full_barcode)

count_matrix    <- count_matrix[, samples_to_keep]
sample_metadata <- sample_metadata[match(samples_to_keep, rownames(sample_metadata)), ]
stopifnot(identical(colnames(count_matrix), rownames(sample_metadata)))

sample_metadata$full_barcode <- colnames(count_matrix)
sample_metadata$patient_id   <- substr(sample_metadata$full_barcode, 1, 12)

preprocess_summary <- bind_rows(preprocess_summary, tibble(
  step = "Deduplicated patients", genes = nrow(count_matrix), samples = ncol(count_matrix)
))
cat("\nPreprocessing so far:\n"); print(preprocess_summary)

# -----------------------------
# Join clinical + PAM50, with EXPLICIT suffixes
# -----------------------------
clinical_patient_col_candidates <- c("submitter_id", "bcr_patient_barcode", "patient_id")
clinical_patient_col <- intersect(clinical_patient_col_candidates, colnames(clinical_brca))[1]
if (is.na(clinical_patient_col)) {
  stop("None of the expected patient-ID columns found in clinical_brca: ",
       paste(colnames(clinical_brca), collapse = ", "))
}
clinical_brca$patient_id <- clinical_brca[[clinical_patient_col]]

pam50_patient_col <- grep("barcode|patient|bcr", colnames(pam50_brca), ignore.case = TRUE, value = TRUE)[1]
pam50_subtype_col <- grep("subtype|pam50",       colnames(pam50_brca), ignore.case = TRUE, value = TRUE)[1]

pam50_clean <- pam50_brca %>%
  mutate(patient_id = substr(.data[[pam50_patient_col]], 1, 12),
         PAM50      = .data[[pam50_subtype_col]]) %>%
  dplyr::select(patient_id, PAM50) %>%
  distinct(patient_id, .keep_all = TRUE)

sample_metadata2 <- sample_metadata %>%
  left_join(clinical_brca, by = "patient_id", suffix = c("_SE", "_GDCclinic")) %>%
  left_join(pam50_clean,   by = "patient_id")

rownames(sample_metadata2) <- sample_metadata2$full_barcode
sample_metadata2 <- sample_metadata2[colnames(count_matrix), ]
stopifnot(identical(rownames(sample_metadata2), colnames(count_matrix)))

cat("\nPAM50 distribution:\n"); print(table(sample_metadata2$PAM50, useNA = "always"))

# Resolve to ONE authoritative, documented clinical column set.
clinical_resolved <- resolve_clinical_columns(sample_metadata2)
sample_metadata2 <- bind_cols(
  sample_metadata2,
  clinical_resolved %>% dplyr::select(-patient_id)
)

cat("\nResolved clinical columns now available on sample_metadata2:\n")
print(setdiff(names(CLINICAL_COLUMN_MAP), "patient_id"))

# -----------------------------
# Integer rounding BEFORE the low-count filter (FIX: previously the
# filter ran on unrounded fractional STAR multi-mapped counts, then
# rounding happened only afterward for DESeq2's benefit -- filtering
# and modeling should see the SAME integer-valued matrix).
# -----------------------------
count_matrix_int <- round(count_matrix)

# -----------------------------
# Low-count gene filter (count >= 10 in >= 10% of samples), using
# ceiling() so the "at least 10%" criterion is never silently relaxed
# by round-half-down behavior on the sample-count threshold.
# -----------------------------
min_samples            <- ceiling(0.10 * ncol(count_matrix_int))
keep_genes              <- rowSums(count_matrix_int >= 10) >= min_samples
count_matrix_filtered   <- count_matrix_int[keep_genes, ]

preprocess_summary <- bind_rows(preprocess_summary, tibble(
  step = "Low-count genes removed",
  genes = nrow(count_matrix_filtered), samples = ncol(count_matrix_filtered)
))
cat("\nFinal preprocessing summary:\n"); print(preprocess_summary)

# -----------------------------
# DESeq2 VST -- UNSUPERVISED (design = ~1, blind = TRUE). Uses ONLY
# the (now-integer) count matrix; clinical/PAM50 columns in colData
# exist for downstream bookkeeping only, NOT as model covariates.
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData   = sample_metadata2,
  design    = ~ 1
)
vst_obj    <- vst(dds, blind = TRUE)
vst_matrix <- assay(vst_obj)
cat("\nVST matrix dimensions:", dim(vst_matrix), "\n")

# -----------------------------
# Top variable genes (ranked AFTER VST)
# -----------------------------
gene_variance <- apply(vst_matrix, 1, var)
top2000_genes <- names(sort(gene_variance, decreasing = TRUE))[1:2000]
top5000_genes <- names(sort(gene_variance, decreasing = TRUE))[1:5000]
vst_top2000   <- vst_matrix[top2000_genes, ]
vst_top5000   <- vst_matrix[top5000_genes, ]

# -----------------------------
# Save outputs
# -----------------------------
saveRDS(count_matrix_filtered, here::here("data/processed/count_matrix_filtered.rds"))
saveRDS(vst_matrix,            here::here("data/processed/vst_matrix_all_genes.rds"))
saveRDS(vst_top2000,           here::here("data/processed/vst_top2000_genes.rds"))
saveRDS(vst_top5000,           here::here("data/processed/vst_top5000_genes.rds"))
saveRDS(sample_metadata2,      here::here("data/processed/sample_metadata_matched.rds"))
saveRDS(gene_annotation,       here::here("data/processed/gene_annotation.rds"))

sample_metadata2_clean <- data.frame(lapply(sample_metadata2, function(x) {
  if (is.list(x)) sapply(x, function(v) paste(unlist(v), collapse = ";")) else x
}), check.names = FALSE)
write_csv(sample_metadata2_clean, here::here("data/processed/sample_metadata_matched.csv"))
write_csv(preprocess_summary,     here::here("results/tables/table_preprocessing_summary.csv"))

write_csv(
  tibble(gene_id = names(gene_variance), variance = gene_variance) %>% arrange(desc(variance)),
  here::here("results/tables/table_gene_variance.csv")
)

writeLines(
  c("NOTE: vst_top2000_genes.rds and vst_top5000_genes.rds are selected by",
    "expression variance WITHIN THIS TCGA-BRCA COHORT ONLY. They are a",
    "TCGA-optimized feature set, not an independently-derived or",
    "literature-curated breast cancer gene signature. External validation",
    "(13_external_validation.R) inherits this TCGA-defined feature space",
    "when projecting METABRIC samples onto TCGA cluster centroids; this",
    "is stated as a limitation in the manuscript."),
  here::here("results/tables/NOTE_tcga_specific_feature_selection.txt")
)

# Persist the unsupervised-discovery design statement explicitly: no
# internal train/holdout split exists anywhere upstream of clustering.
writeLines(
  c("DESIGN NOTE: unsupervised discovery, no internal holdout split.",
    "================================================================",
    "Low-count gene filtering, top-variance feature selection, SNF",
    "clustering, and cluster-characterization analyses (DE, pathway,",
    "master-regulator) all draw on the FULL TCGA-BRCA discovery cohort",
    "with no train/test partition. This is standard and defensible for",
    "an unsupervised molecular subtyping study (there is no outcome",
    "label being predicted that would require a holdout), but it means",
    "cluster-characterization p-values (scripts 08, 09, 11, 12) reflect",
    "testing AFTER unsupervised grouping on the SAME data used to form",
    "the groups, and should be interpreted as descriptive/hypothesis-",
    "generating within TCGA, with independent replication carried by",
    "the METABRIC external validation (script 13) alone."),
  here::here("results/tables/NOTE_unsupervised_discovery_design.txt")
)

log_session_info(script_name, key_packages = c("DESeq2", "SummarizedExperiment"))
cat("\n\u2713 03_preprocess_expression.R complete.\n")

close_logger(log_con, script_name)
