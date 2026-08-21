# ============================================================
# 12_master_regulator_network.R  (publication-ready revision, FINAL)
# Cluster-specific master regulator TFs + simplified TF->target network
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(dorothea)
  library(igraph)
  library(ggraph)
  library(pheatmap)
  library(purrr)
})
source(here::here("R", "utils.R"))

script_name <- "12_master_regulator_network"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

circularity_note <- c(
  "METHODOLOGICAL CAVEAT: TF 'master regulator' identification is partly",
  "circular by construction",
  "========================================================================",
  "TF activity (VIPER/DoRothEA) is one of the TWO feature views fused by SNF",
  "to PRODUCE the cluster labels tested in this script. Finding that TF",
  "activity differs significantly across SNF clusters is therefore expected",
  "to a substantial degree simply because SNF was designed to find groups",
  "that separate in (among other things) TF-activity space -- it is not, on",
  "its own, independent evidence of true master-regulator biology.",
  "",
  "INTERPRETATION: results here should be corroborated by (a) standardized",
  "effect size (a genuine one-vs-rest Cohen's d, not just statistical",
  "significance, which is near-guaranteed for SOME TFs by construction),",
  "(b) concordance with independently-derived marker-gene/pathway results",
  "(script 11), and (c) consistency with established breast cancer subtype",
  "biology in the literature (e.g. ESR1/FOXA1/GATA3 activity in luminal-like",
  "clusters), rather than treated as a novel statistical discovery in isolation."
)
writeLines(circularity_note, here::here("results/tables/NOTE_tf_circularity_caveat.txt"))
cat(circularity_note, sep = "\n"); cat("\n")

tf_activity      <- readRDS(here::here("data/processed/tf_activity_viper_AC_primary.rds"))
final_labels_obj <- readRDS(here::here("results/objects/final_snf_cluster_labels.rds"))
de_table         <- read_csv(here::here("results/tables/table_de_markers_k2_k3.csv"), show_col_types = FALSE)

data(dorothea_hs, package = "dorothea")
regulon_ac <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))

run_master_regulators_for_k <- function(cluster_labels_k, k) {
  cat("\n========== Master regulator analysis at k =", k, "==========\n")
  
  aligned <- align_samples(tf = tf_activity,
                           meta = cluster_labels_k %>% column_to_rownames("full_barcode"))
  tf_act_k <- aligned$tf
  meta_k   <- aligned$meta
  
  tf_long <- as.data.frame(t(tf_act_k)) %>%
    rownames_to_column("full_barcode") %>%
    pivot_longer(-full_barcode, names_to = "TF", values_to = "activity") %>%
    left_join(meta_k %>% rownames_to_column("full_barcode") %>% dplyr::select(full_barcode, SNF_cluster),
              by = "full_barcode")
  
  # --- Kruskal-Wallis per TF -- FIX: single kruskal.test() call per TF,
  # both kw_stat and kw_p extracted from the same result object. ---
  
  tf_test_table <- tf_long %>%
    group_by(TF) %>%
    summarise(
      kw_result = list(tryCatch(kruskal.test(activity ~ SNF_cluster), error = function(e) NULL)),
      .groups = "drop"
    ) %>%
    mutate(
      kw_stat = map_dbl(kw_result, ~ if (is.null(.x)) NA_real_ else .x$statistic),
      kw_p    = map_dbl(kw_result, ~ if (is.null(.x)) NA_real_ else .x$p.value)
    ) %>%
    dplyr::select(-kw_result) %>%
    mutate(FDR = p.adjust(kw_p, method = "BH"), k = k)
  # --- GENUINE one-vs-rest Cohen's d per (TF, cluster), pooled SD ---
  cluster_levels <- sort(unique(tf_long$SNF_cluster))
  tf_names <- unique(tf_long$TF)
  
  effect_size_grid <- expand_grid(SNF_cluster = cluster_levels, TF = tf_names)
  
  tf_wide_by_tf <- split(tf_long, tf_long$TF)
  cohens_d_table <- map_dfr(tf_names, function(tfx) {
    d <- tf_wide_by_tf[[tfx]]
    map_dfr(cluster_levels, function(cl) {
      tibble(
        TF = tfx, SNF_cluster = cl,
        mean_in = mean(d$activity[d$SNF_cluster == cl], na.rm = TRUE),
        mean_rest = mean(d$activity[d$SNF_cluster != cl], na.rm = TRUE),
        cohens_d = cohens_d_one_vs_rest(d$activity, d$SNF_cluster, cl)
      )
    })
  })
  
  tf_cluster_summary <- cohens_d_table %>%
    mutate(activity_difference = mean_in - mean_rest) %>%
    left_join(tf_test_table, by = "TF") %>%
    filter(FDR < 0.05) %>%
    arrange(SNF_cluster, desc(abs(cohens_d)))
  
  top10_per_cluster <- tf_cluster_summary %>%
    group_by(SNF_cluster) %>%
    slice_max(order_by = abs(cohens_d), n = 10) %>%
    ungroup() %>%
    mutate(k = k)
  
  cat("Top master regulators (by GENUINE one-vs-rest Cohen's d) per cluster:\n")
  print(top10_per_cluster %>% dplyr::select(SNF_cluster, TF, activity_difference, cohens_d, FDR))
  
  # --- Triangulation with script 11 DE results ---
  # FIX: join is now keyed on BOTH gene symbol AND cluster, so each
  # master regulator maps to exactly one DE row -- the result from the
  # SAME cluster it was identified as a master regulator for -- rather
  # than fanning out across every cluster tested for that gene.
  # NOTE: rename() happens BEFORE the join here, so SYMBOL no longer
  # exists as a column name in this table by the time left_join() runs.
  # The join key must therefore reference "TF" (not "SYMBOL") on both
  # sides for the gene-symbol match.
  de_k <- de_table %>% filter(k == !!k)
  triangulation <- top10_per_cluster %>%
    left_join(
      de_k %>% dplyr::select(SYMBOL, cluster, log2FoldChange, padj, used_in_clustering) %>%
        rename(TF = SYMBOL, de_log2FC = log2FoldChange, de_padj = padj),
      by = c("TF", "SNF_cluster" = "cluster")
    ) %>%
    mutate(de_concordant_cluster = !is.na(de_log2FC))
  
  cat("\nTriangulation with script 11 DE results (is the top TF gene itself\n",
      "differentially expressed in the SAME cluster it was identified as a\n",
      "master regulator for? `used_in_clustering` flags whether that gene was\n",
      "also an SNF input feature, i.e. how independent this corroboration is):\n", sep = "")
  print(triangulation %>% dplyr::select(SNF_cluster, TF, cohens_d, de_concordant_cluster, de_log2FC, de_padj, used_in_clustering))
  
  # Sanity check: this join must NOT fan out. Each top10 row should map
  # to at most one DE row (same-cluster match), so triangulation should
  # have the SAME number of rows as top10_per_cluster.
  if (nrow(triangulation) != nrow(top10_per_cluster)) {
    warning("k=", k, ": triangulation row count (", nrow(triangulation),
            ") does not match top10_per_cluster row count (", nrow(top10_per_cluster),
            ") -- the join may be fanning out again. Investigate before trusting this table.")
  } else {
    cat("\nJoin integrity check passed: triangulation has exactly", nrow(triangulation),
        "rows, matching top10_per_cluster (no fanout).\n")
  }
  
  list(top10 = top10_per_cluster, triangulation = triangulation,
       tf_test_table = tf_test_table, tf_act_k = tf_act_k, meta_k = meta_k)
  }

res_k2 <- run_master_regulators_for_k(final_labels_obj$k2, 2)
res_k3 <- run_master_regulators_for_k(final_labels_obj$k3, 3)

top10_combined <- bind_rows(res_k2$top10, res_k3$top10)
triangulation_combined <- bind_rows(res_k2$triangulation, res_k3$triangulation)

write_csv(top10_combined,         here::here("results/tables/table_master_regulators_k2_k3.csv"))
write_csv(triangulation_combined, here::here("results/tables/table_master_regulator_triangulation_k2_k3.csv"))

# -----------------------------
# Network visualization (literature-curated DoRothEA edges only)
# -----------------------------
plot_network_for_cluster <- function(top10_df, regulon, k, cl) {
  tfs_here <- top10_df %>% filter(SNF_cluster == cl) %>% pull(TF)
  edges <- regulon %>%
    filter(tf %in% tfs_here) %>%
    # Break ties in |mor| by DoRothEA confidence tier (A > B > C, i.e.
    # alphabetical order = quality order here), then by target name for
    # full reproducibility -- rather than relying on with_ties=FALSE's
    # arbitrary original-row-order tie-break.
    arrange(tf, confidence, desc(abs(mor)), target) %>%
    group_by(tf) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    transmute(from = tf, to = target, mor = mor)

  
  if (nrow(edges) == 0) return(invisible(NULL))
  g <- igraph::graph_from_data_frame(edges, directed = TRUE)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(color = mor), arrow = arrow(length = unit(2, "mm")), alpha = 0.6) +
    scale_edge_color_gradient2(low = "blue", mid = "grey80", high = "red", midpoint = 0) +
    geom_node_point(size = 3) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    labs(title = paste0("Top regulators, cluster ", cl, " (k=", k,
                        ") -- literature-curated DoRothEA edges, top-5 |mor| targets per TF")) +
    theme_void()
  ggsave(here::here("results/figures",
                    paste0("figure_tf_network_k", k, "_", cl, ".png")),
         p, width = 7, height = 6, dpi = 300)
}

for (cl in unique(res_k3$top10$SNF_cluster)) plot_network_for_cluster(res_k3$top10, regulon_ac, 3, cl)
for (cl in unique(res_k2$top10$SNF_cluster)) plot_network_for_cluster(res_k2$top10, regulon_ac, 2, cl)

log_session_info(script_name, key_packages = c("dorothea", "igraph", "ggraph"))
cat("\n\u2713 12_master_regulator_network.R complete. Genuine one-vs-rest Cohen's d used for ranking;",
    "triangulation join fixed to match on cluster (no fanout); Kruskal-Wallis computed once per TF.\n")

close_logger(log_con, script_name)