# ============================================================
# 00_setup.R  (publication-ready revision)
# Project setup for: TF Activity-Integrated SNF for BRCA Subtyping
# ============================================================
# CHANGES in this revision (audit follow-up):
#   - Critical, API-sensitive packages (decoupleR, dorothea, msigdbr,
#     TCGAbiolinks) now STOP the pipeline if below the pinned minimum
#     version, rather than only warning. A version mismatch in these
#     packages can silently change NUMERICAL results (regulon content,
#     score semantics), not just column names -- a warning is not a
#     strong enough guardrail for that class of risk.
#   - renv snapshot guidance corrected: a single renv::snapshot() call
#     immediately after this script only captures packages loaded by
#     00_setup.R itself. This revision instructs (and provides a
#     helper) to snapshot AFTER all 13 scripts have run at least once,
#     using renv::dependencies() scanned across scripts/, so the
#     lockfile captures every package actually used by the pipeline
#     (survRM2, aricode, ggraph, ReactomePA, etc.), not just setup-time
#     packages.
# ============================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# ---- 0. Pin the Bioconductor release used for this project -------
BIOC_VERSION <- "3.19"

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = BIOC_VERSION, ask = FALSE, update = FALSE)

# ---- 1. CRAN packages ----------------------------------------------
cran_packages <- c(
  "here", "tidyverse", "ggplot2", "pheatmap", "RColorBrewer",
  "cluster", "survival", "survminer", "igraph", "ggraph",
  "mclust", "aricode", "renv", "broom", "survRM2"
)
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

# ---- 2. Bioconductor packages --------------------------------------
bioc_packages <- c(
  "TCGAbiolinks", "SummarizedExperiment", "DESeq2",
  "AnnotationDbi", "org.Hs.eg.db", "dorothea", "decoupleR",
  "limma", "clusterProfiler", "ReactomePA", "msigdbr",
  "enrichplot"
)
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# ---- 3. SNFtool (not always on the current CRAN mirror) ------------
if (!requireNamespace("SNFtool", quietly = TRUE)) {
  tryCatch(
    install.packages("SNFtool"),
    error = function(e) {
      message("SNFtool not found on CRAN mirror; trying CRAN archive...")
      install.packages(
        "https://cran.r-project.org/src/contrib/Archive/SNFtool/SNFtool_2.3.1.tar.gz",
        repos = NULL, type = "source"
      )
    }
  )
}

# ---- 4. HARD minimum-version checks for API-sensitive packages -----
# These packages have had breaking changes to ARGUMENT NAMES, OUTPUT
# COLUMNS, or (critically) SCORE SEMANTICS in the recent past. A
# version mismatch here can silently change numbers, not just crash --
# so this is a stop(), not a warning().
check_min_version <- function(pkg, min_version, hard_stop = TRUE) {
  inst <- tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA)
  if (is.na(inst)) {
    msg <- paste0(pkg, " is not installed; cannot verify minimum version ", min_version)
    if (hard_stop) stop(msg) else { warning(msg); return(invisible(FALSE)) }
  }
  ok <- utils::compareVersion(inst, min_version) >= 0
  if (!ok) {
    msg <- paste0(
      pkg, " version ", inst, " is older than the minimum tested/pinned version ",
      min_version, ". This package's output schema and/or score semantics have ",
      "changed across versions historically (regulon format, VIPER score columns, ",
      "msigdbr category->collection rename). Install the pinned version before proceeding."
    )
    if (hard_stop) stop(msg) else warning(msg)
  } else {
    cat("OK:", pkg, inst, ">=", min_version, "\n")
  }
  invisible(ok)
}

check_min_version("msigdbr",      "10.0.0", hard_stop = TRUE)
check_min_version("decoupleR",    "2.4.0",  hard_stop = TRUE)
check_min_version("TCGAbiolinks", "2.31.0", hard_stop = TRUE)
check_min_version("dorothea",     "1.14.0", hard_stop = TRUE)

# ---- 5. Project folder scaffold ------------------------------------
folders <- c(
  "data/raw", "data/processed", "data/external",
  "scripts", "R", "results/tables", "results/figures",
  "results/objects", "logs", "logs/sessioninfo", "thesis"
)
for (f in folders) dir.create(here::here(f), recursive = TRUE, showWarnings = FALSE)

cat(
  "\nNOTE: the following external files must be placed manually before\n",
  "running the corresponding scripts (they are not auto-downloaded):\n",
  "  data/external/tumor_purity.csv         (used by 04_quality_control.R)\n",
  "  data/external/metabric_expression.txt  (required by 13_external_validation.R)\n",
  "  data/external/metabric_clinical.txt    (required by 13_external_validation.R)\n\n",
  sep = ""
)

# ---- 6. Save sessionInfo() and initialize renv ----------------------
writeLines(capture.output(sessionInfo()), here::here("logs", "sessioninfo", "00_setup_sessionInfo.txt"))

if (!file.exists(here::here("renv.lock"))) {
  message("No renv.lock found -- initializing renv for this project now.")
  renv::init(bare = TRUE, restart = FALSE)
}

#' Full-pipeline dependency snapshot. Call this ONCE, AFTER running all
#' 13 scripts at least one time each (e.g. from a final "run_all.R"
#' driver), NOT immediately after 00_setup.R alone -- a snapshot taken
#' here would only capture 00_setup.R's own library() calls and could
#' under-represent packages first loaded by later scripts (survRM2,
#' aricode, ggraph, ReactomePA, msigdbr, etc.), which would make the
#' lockfile an incomplete reproducibility artifact.
snapshot_full_pipeline <- function() {
  deps <- tryCatch(
    renv::dependencies(here::here("scripts"))$Package,
    error = function(e) NULL
  )
  renv::snapshot(prompt = FALSE, packages = deps)
}

cat(
  "\nNOTE: run snapshot_full_pipeline() (defined in this script) AFTER all\n",
  "13 numbered scripts have been executed at least once, so renv.lock\n",
  "captures every package actually used by the full pipeline.\n"
)

cat("\n\u2713 00_setup.R complete. Environment ready.\n")
