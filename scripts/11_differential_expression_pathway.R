# ============================================================
# 11_differential_expression_pathway.R  (publication-ready revision, FINAL)
# Cluster-specific marker genes (DESeq2, one-vs-rest) + GSEA
# (Hallmark, Reactome)
# ============================================================
# CHANGES in this revision (audit follow-up):
#   - map_to_symbol_dedup() is now called with an EXPLICIT
#     `variance_reference = vst_matrix_all_genes`, the SAME reference
#     used in script 05 -- previously this script tie-broke duplicate
#     Ensembl IDs by variance of RAW COUNTS while script 05 tie-broke
#     by variance of VST values, which could silently select a
#     DIFFERENT representative Ensembl ID for the same gene symbol in
#     the two scripts despite both claiming "the same shared rule."
#   - Every marker-gene row is now flagged `used_in_clustering`
#     (TRUE if that gene's symbol was one of the top-2000-variance
#     genes fed into SNF). This directly addresses the audit finding
#     that the persisted circularity caveat treated all DE results as
#     equally circular, when in fact genes that were NOT clustering
#     inputs provide substantially stronger (less mechanically
#     circular) evidence than genes that were. The caveat text is
#     updated to make this distinction explicit.
#   - FIX (this revision): GSEA/gsePathway ranked lists are now built
#     from the FULL DESeq2 result (res_df_full), not from the
#     padj-filtered marker table (res_df). DESeq2's independent
#     filtering sets padj = NA for low-mean-count genes but usually
#     leaves `stat` intact for those genes; the previous version
#     silently dropped every independent-filtering-excluded gene from
#     the pre-ranked GSEA input, which is a deviation from standard
#     pre-ranked GSEA practice (rank on the full tested set) and can
#     shift enrichment results beyond just "fewer genes considered."
#     The padj-filtered table is still used, unchanged, for the
#     marker-gene output table and volcano plots.
#   - FIX (this revision): Reactome ENTREZID de-duplication now tie-
#     breaks by descending |stat| before distinct(), so the retained
#     symbol per Entrez ID is the more extreme one on the SAME
#     statistic used for ranking -- previously the implicit tie-break
#     was "whichever symbol had smallest padj" (an artifact of
#     res_df's row order), which is inconsistent with ranking by stat.
#   - FIX (this revision): results() now called with alpha = 0.05 to
#     match this script's actual significance threshold (padj < 0.05),
#     so DESeq2's independent-filtering threshold is calibrated to the
#     cutoff actually used, rather than the default alpha = 0.1.
#   - FIX (this revision): volcano plot -log10(padj) is clipped at 50
#     so genes with padj rounding to exactly 0 (common at this sample
#     size) are still plotted rather than silently dropped as Inf.
#   - NOTE (this revision): GSEA()/gsePathway() are called with
#     seed = TRUE. In several clusterProfiler versions this fixes an
#     INTERNAL seed for the permutation step rather than deferring to
#     the caller's RNG state -- i.e. it may NOT be governed by the
#     set_pipeline_seed(offset = 1100) call below. This does not
#     affect within-script reproducibility (results are still
#     deterministic run-to-run), but the offset should not be assumed
#     to control GSEA permutation randomness specifically. Verify
#     against your installed clusterProfiler version's documentation
#     if this matters for your methods write-up.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(ReactomePA)
  library(msigdbr)
})
source(here::here("R", "utils.R"))

script_name <- "11_differential_expression_pathway"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

count_matrix_filtered <- readRDS(here::here("data/processed/count_matrix_filtered.rds"))
vst_matrix_all_genes  <- readRDS(here::here("data/processed/vst_matrix_all_genes.rds"))
vst_top2000           <- readRDS(here::here("data/processed/vst_top2000_genes.rds"))
sample_metadata       <- readRDS(here::here("data/processed/sample_metadata_matched.rds"))
final_labels_obj      <- readRDS(here::here("results/objects/final_snf_cluster_labels.rds"))

circularity_note <- c(
  "METHODOLOGICAL CAVEAT: differential expression after unsupervised clustering",
  "================================================================================",
  "Marker-gene (DESeq2) and pathway-enrichment (GSEA) p-values in this script are",
  "computed using the SAME expression data used to derive the SNF cluster labels",
  "being tested ('testing after clustering'). DESeq2's p-values assume the",
  "compared groups are fixed independently of the tested data; that assumption",
  "is violated here.",
  "",
  "IMPORTANT DISTINCTION (added in this revision): circularity is NOT uniform",
  "across all reported markers. Genes that were themselves INPUT FEATURES to",
  "SNF clustering (the top-2000-variance gene set) are expected to differ",
  "across clusters substantially BY CONSTRUCTION. Genes NOT used as clustering",
  "features (the majority of genes tested here) differing across clusters is a",
  "DOWNSTREAM CONSEQUENCE of the clustering, not a direct mechanical artifact --",
  "this is meaningfully stronger, less circular evidence. Every row in",
  "table_de_markers_k2_k3.csv carries a `used_in_clustering` flag so a reader",
  "can distinguish these two evidence tiers directly. Fold-change rankings and",
  "enrichment DIRECTION remain informative for BOTH tiers; raw/padj p-value",
  "magnitudes should not be over-interpreted for genes with used_in_clustering",
  "= TRUE. Partial mitigation comes from external replication (METABRIC,",
  "script 13), though gene-level DE/pathway results are not independently",
  "re-tested in that external cohort.",
  "",
  "GSEA RANKING NOTE (added in this revision): pre-ranked GSEA (Hallmark and",
  "Reactome) is run on the FULL tested gene set (ranked by DESeq2's `stat`),",
  "not restricted to genes surviving padj-based independent filtering. This",
  "follows standard pre-ranked GSEA practice and avoids truncating the ranked",
  "list before it reaches the enrichment step."
)
writeLines(circularity_note, here::here("results/tables/NOTE_de_after_clustering_caveat.txt"))
cat(circularity_note, sep = "\n"); cat("\n")

# -----------------------------
# Symbol mapping -- SAME variance reference (VST matrix) as script 05,
# regardless of the fact that this script's returned matrix is raw
# counts. This is the fix: the tie-break rule is now scale-consistent
# across the whole pipeline.
# -----------------------------
mapped <- map_to_symbol_dedup(
  count_matrix_filtered,
  variance_reference = vst_matrix_all_genes,
  org_db = org.Hs.eg.db
)
count_symbol_mat <- round(mapped$matrix)
symbol_map       <- mapped$map

# FIX: rownames(vst_top2000) are Ensembl IDs (set in script 03, before
# symbol mapping ever occurs), but res_df$SYMBOL is a gene SYMBOL --
# comparing them directly always returns FALSE. Map vst_top2000's
# Ensembl rownames to symbols using the SAME map_to_symbol_dedup()
# already computed above (`symbol_map`), so the clustering-input
# comparison uses a consistent gene identity namespace.
clustering_input_ensembl <- rownames(vst_top2000)
clustering_input_symbols <- symbol_map %>%
  filter(ENSEMBL %in% sub("\\..*$", "", clustering_input_ensembl)) %>%
  pull(SYMBOL) %>%
  unique()

cat("Clustering-input gene count (Ensembl):", length(clustering_input_ensembl),
    " -> resolved to", length(clustering_input_symbols), "unique symbols\n",
    "(some loss expected: this script's symbol_map was built from the FULL",
    "filtered gene set, dedup'd by variance -- a top-2000-variance Ensembl ID",
    "could map to a symbol whose HIGHEST-variance representative elsewhere",
    "was a DIFFERENT Ensembl ID, in which case it won't appear as a symbol here.",
    "This is a known, minor, and correctly-behaved edge case of the dedup rule,",
    "not a bug.)\n")


hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

run_de_and_gsea_for_k <- function(cluster_labels_k, k) {
  cat("\n========== DE + GSEA at k =", k, "==========\n")
  
  meta_k <- sample_metadata %>%
    inner_join(cluster_labels_k, by = "full_barcode") %>%
    column_to_rownames("full_barcode")
  meta_k <- meta_k[colnames(count_symbol_mat), , drop = FALSE]
  stopifnot(identical(rownames(meta_k), colnames(count_symbol_mat)))
  
  cluster_levels <- sort(unique(meta_k$SNF_cluster))
  
  de_results_all <- list(); gsea_hallmark_all <- list(); gsea_reactome_all <- list()
  
  for (cl in cluster_levels) {
    cat("\n--- Cluster", cl, "(k=", k, ") vs rest ---\n")
    
    meta_k$cluster_vs_rest <- factor(ifelse(meta_k$SNF_cluster == cl, "InCluster", "Other"),
                                     levels = c("Other", "InCluster"))
    
    dds <- DESeqDataSetFromMatrix(
      countData = count_symbol_mat, colData = meta_k, design = ~ cluster_vs_rest
    )
    dds <- DESeq(dds, quiet = TRUE)
    # FIX: alpha = 0.05 matches this script's actual significance cutoff
    # (padj < 0.05 below), rather than DESeq2's default alpha = 0.1,
    # so independent filtering is calibrated to the threshold we use.
    res <- results(dds, contrast = c("cluster_vs_rest", "InCluster", "Other"), alpha = 0.05)
    
    # FULL result set (no padj filtering) -- this is the correct input
    # for pre-ranked GSEA. Kept separate from the padj-filtered marker
    # table below.
    res_df_full <- as.data.frame(res) %>%
      rownames_to_column("SYMBOL") %>%
      mutate(used_in_clustering = SYMBOL %in% clustering_input_symbols)
    
    # padj-filtered table -- used for the marker-gene output table and
    # volcano plot ONLY. Not used for GSEA ranking (see fix above).
    res_df <- res_df_full %>%
      filter(!is.na(padj)) %>%
      arrange(padj)
    
    res_df$k <- k; res_df$cluster <- cl
    de_results_all[[cl]] <- res_df
    
    sig_markers <- res_df %>% filter(padj < 0.05, abs(log2FoldChange) >= 1)
    cat("Significant markers (padj<0.05, |LFC|>=1):", nrow(sig_markers),
        " (of which", sum(sig_markers$used_in_clustering), "were clustering-input genes)\n")
    cat("Genes excluded from padj (independent filtering):",
        sum(is.na(res_df_full$padj)),
        "-- these ARE still included in the GSEA ranked list via res_df_full\n")
    
    # FIX: clip -log10(padj) at 50 so genes with padj rounding to exactly
    # 0 (common at this sample size) still plot instead of silently
    # dropping as Inf.
    p_volcano <- ggplot(res_df, aes(log2FoldChange, pmin(-log10(padj), 50), shape = used_in_clustering)) +
      geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) >= 1), alpha = 0.5, size = 0.8) +
      scale_color_manual(values = c("grey70", "firebrick"), guide = "none") +
      labs(title = paste0("Cluster ", cl, " vs rest (k=", k, ")"),
           subtitle = "Shape: triangle = gene WAS an SNF clustering input feature",
           x = "log2 fold change", y = "-log10(padj)  [clipped at 50]") +
      theme_minimal()
    ggsave(here::here("results/figures",
                      paste0("figure_volcano_k", k, "_", cl, ".png")),
           p_volcano, width = 6, height = 5, dpi = 300)
    
    # FIX: ranked list now built from res_df_full (the FULL tested set,
    # not padj-filtered), consistent with standard pre-ranked GSEA
    # practice.
    ranked_df <- res_df_full %>%
      filter(!is.na(stat)) %>%
      arrange(desc(stat))
    ranked <- ranked_df$stat
    names(ranked) <- ranked_df$SYMBOL
    cat("GSEA (Hallmark) ranked list size:", length(ranked), "genes\n")
    
    gsea_h <- tryCatch(
      GSEA(ranked, TERM2GENE = hallmark_sets, pvalueCutoff = 1, seed = TRUE, verbose = FALSE),
      error = function(e) { message("Hallmark GSEA failed for cluster ", cl, ": ", conditionMessage(e)); NULL }
    )
    if (!is.null(gsea_h)) {
      gh_df <- as.data.frame(gsea_h) %>% mutate(k = k, cluster = cl)
      gsea_hallmark_all[[cl]] <- gh_df
    }
    
    # FIX: also built from res_df_full now. Tie-break for duplicate
    # ENTREZIDs is by descending |stat| (the statistic actually used
    # for ranking), not by whatever row order res_df happened to have
    # from its padj sort.
    entrez_ranked_df <- res_df_full %>%
      inner_join(symbol_map, by = "SYMBOL") %>%
      filter(!is.na(ENTREZID), !is.na(stat)) %>%
      arrange(desc(abs(stat))) %>%
      distinct(ENTREZID, .keep_all = TRUE) %>%
      arrange(desc(stat))
    ranked_entrez <- entrez_ranked_df$stat
    names(ranked_entrez) <- entrez_ranked_df$ENTREZID
    cat("GSEA (Reactome) ranked list size:", length(ranked_entrez), "genes\n")
    
    gsea_r <- tryCatch(
      gsePathway(ranked_entrez, organism = "human", pvalueCutoff = 1, seed = TRUE, verbose = FALSE),
      error = function(e) { message("Reactome GSEA failed for cluster ", cl, ": ", conditionMessage(e)); NULL }
    )
    if (!is.null(gsea_r)) {
      gr_df <- as.data.frame(gsea_r) %>% mutate(k = k, cluster = cl)
      gsea_reactome_all[[cl]] <- gr_df
    }
  }
  
  list(
    de = bind_rows(de_results_all),
    gsea_hallmark = bind_rows(gsea_hallmark_all),
    gsea_reactome = bind_rows(gsea_reactome_all)
  )
}

set_pipeline_seed(offset = 1100)
res_k2 <- run_de_and_gsea_for_k(final_labels_obj$k2, 2)
set_pipeline_seed(offset = 1100)
res_k3 <- run_de_and_gsea_for_k(final_labels_obj$k3, 3)

de_table        <- bind_rows(res_k2$de, res_k3$de)
gsea_hallmark   <- bind_rows(res_k2$gsea_hallmark, res_k3$gsea_hallmark)
gsea_reactome   <- bind_rows(res_k2$gsea_reactome, res_k3$gsea_reactome)

write_csv(de_table,      here::here("results/tables/table_de_markers_k2_k3.csv"))
write_csv(gsea_hallmark, here::here("results/tables/table_gsea_hallmark_k2_k3.csv"))
write_csv(gsea_reactome, here::here("results/tables/table_gsea_reactome_k2_k3.csv"))

log_session_info(script_name, key_packages = c("DESeq2", "clusterProfiler", "ReactomePA", "msigdbr"))
cat("\n\u2713 11_differential_expression_pathway.R complete. used_in_clustering flag added",
    "for circularity transparency; GSEA now ranks on the full tested gene set.\n")

close_logger(log_con, script_name)