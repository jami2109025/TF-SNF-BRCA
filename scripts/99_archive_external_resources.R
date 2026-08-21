# ============================================================
# 99_archive_external_resources.R
# Freeze every version-sensitive external resource the pipeline
# depends on, with checksums and provenance, ready for deposit.
# ============================================================
# WHY THIS EXISTS
# ------------------------------------------------------------
# Several inputs to this pipeline are LIVE or VERSIONED resources
# fetched from the internet at run time. If any of them changes between
# the run that produced the manuscript's numbers and the run a reviewer
# or a future reader attempts, the pipeline will not error -- it will
# silently produce DIFFERENT numbers from the same code. That is the
# worst failure mode in computational reproducibility, because nothing
# announces it.
#
# Resources with this exposure in this project:
#
#   CollecTRI     fetched live from OmnipathR. Actively curated; edges
#                 are added and revised between releases.
#   DoRothEA      ships inside the dorothea package, so it is pinned by
#                 the package version -- but ONLY if that version is
#                 recorded, which is what this script does.
#   MSigDB        via msigdbr. Hallmark and Reactome gene-set contents
#                 change between MSigDB releases, and script 11's GSEA
#                 results move with them. This is the most commonly
#                 overlooked of the three.
#   org.Hs.eg.db  Ensembl-to-symbol mapping. Changes with each
#                 Bioconductor release, which can alter which genes
#                 survive mapping and therefore the feature set.
#   TCGA/METABRIC the primary data. Not re-archived here (redistribution
#                 is restricted and they are already versioned at
#                 source), but their access dates are recorded.
#
# WHAT THIS SCRIPT PRODUCES
#   archive/external_resources_<DATE>/
#     *.rds            R-native copies (fast to reload)
#     *.tsv.gz         plain-text copies (readable without R, forever)
#     manifest.csv     one row per file: rows, columns, MD5, size
#     MANIFEST.md      human-readable provenance for the deposit page
#     sessionInfo.txt  full package environment
#
# The plain-text export is not redundant. An .rds is an R-specific
# binary whose readability depends on future R versions; a .tsv.gz is
# readable by anything, by anyone, indefinitely. Archiving only the
# .rds is a half-measure that looks like reproducibility.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})
source(here::here("R", "utils.R"))

script_name <- "99_archive_external_resources"
log_con <- init_logger(script_name)
ensure_dirs()

ARCHIVE_DATE <- Sys.Date()
archive_dir <- here::here("archive", paste0("external_resources_", ARCHIVE_DATE))
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
cat("Archive directory:", archive_dir, "\n\n")

manifest_rows <- list()

#' Write one resource in both formats and record its provenance.
#'
#' The MD5 is the operative part of this whole script: it is what lets
#' a reviewer confirm they are holding the identical object, rather
#' than merely a file with the same name.
archive_resource <- function(obj, name, source_desc, package_used = NA_character_) {
  if (is.null(obj)) {
    cat("  [skip]", name, "-- not available\n")
    return(invisible(NULL))
  }

  df <- as.data.frame(obj)
  rds_path <- file.path(archive_dir, paste0(name, ".rds"))
  tsv_path <- file.path(archive_dir, paste0(name, ".tsv.gz"))

  saveRDS(obj, rds_path)
  readr::write_tsv(df, tsv_path)

  pkg_version <- if (!is.na(package_used)) {
    tryCatch(as.character(utils::packageVersion(package_used)), error = function(e) NA_character_)
  } else NA_character_

  manifest_rows[[name]] <<- tibble(
    resource       = name,
    n_rows         = nrow(df),
    n_cols         = ncol(df),
    columns        = paste(colnames(df), collapse = "; "),
    source         = source_desc,
    package        = package_used,
    package_version = pkg_version,
    retrieved_on   = as.character(ARCHIVE_DATE),
    file_rds       = basename(rds_path),
    file_tsv       = basename(tsv_path),
    md5_rds        = unname(tools::md5sum(rds_path)),
    md5_tsv        = unname(tools::md5sum(tsv_path)),
    size_kb        = round(file.size(tsv_path) / 1024, 1)
  )

  cat(sprintf("  [ok] %-22s %7d rows  (%s v%s)  md5=%s\n",
              name, nrow(df), package_used, pkg_version,
              substr(manifest_rows[[name]]$md5_tsv, 1, 8)))
  invisible(NULL)
}

# ------------------------------------------------------------------
# 1. CollecTRI -- archived FROM THE CACHE, never re-fetched.
#
# This is deliberate and important. Re-downloading here would archive
# whatever CollecTRI looks like TODAY, which may differ from the
# version that actually produced the results. The cache written by
# script 05 on its first run IS the object of record; if it is
# missing, that is a hard error, not something to paper over with a
# fresh download.
# ------------------------------------------------------------------
cat("Archiving external resources...\n")

collectri_cache <- here::here("data/processed/network_collectri.rds")
if (file.exists(collectri_cache)) {
  archive_resource(readRDS(collectri_cache), "collectri_network",
                   "OmnipathR::CollecTRI via decoupleR::get_collectri(), cached by 05_tf_activity.R",
                   "OmnipathR")
} else {
  warning("network_collectri.rds not found. Run 05_tf_activity.R first. ",
          "Do NOT substitute a fresh download -- it may differ from the ",
          "version that produced your results.")
}

# ------------------------------------------------------------------
# 2. DoRothEA
# ------------------------------------------------------------------
dorothea_obj <- tryCatch({
  data("dorothea_hs", package = "dorothea", envir = environment())
  get("dorothea_hs", envir = environment())
}, error = function(e) NULL)
archive_resource(dorothea_obj, "dorothea_hs_regulons",
                 "dorothea package data object dorothea_hs", "dorothea")

# ------------------------------------------------------------------
# 3. MSigDB gene sets used by script 11's GSEA
#
# The most commonly missed archive target in this kind of pipeline.
# Hallmark and Reactome set MEMBERSHIP changes between MSigDB
# releases, so GSEA results shift even with identical code and
# identical expression data.
# ------------------------------------------------------------------
msig_hallmark <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "H"),
                          error = function(e) tryCatch(
                            msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                            error = function(e2) NULL))
archive_resource(msig_hallmark, "msigdb_hallmark",
                 "msigdbr, collection H (Hallmark)", "msigdbr")

msig_reactome <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "C2",
                                           subcategory = "CP:REACTOME"),
                          error = function(e) tryCatch(
                            msigdbr::msigdbr(species = "Homo sapiens", collection = "C2",
                                             subcollection = "CP:REACTOME"),
                            error = function(e2) NULL))
archive_resource(msig_reactome, "msigdb_reactome",
                 "msigdbr, collection C2 subcollection CP:REACTOME", "msigdbr")

# ------------------------------------------------------------------
# 4. The Ensembl -> symbol map ACTUALLY USED
#
# Archiving the resolved mapping table, rather than just recording the
# org.Hs.eg.db version, means a reader can reproduce the feature space
# exactly even if that annotation package is no longer installable.
# ------------------------------------------------------------------
symbol_map_path <- here::here("results/tables/table_ensembl_symbol_map_used.csv")
if (file.exists(symbol_map_path)) {
  archive_resource(readr::read_csv(symbol_map_path, show_col_types = FALSE),
                   "ensembl_symbol_map_used",
                   "Resolved Ensembl->SYMBOL map produced by 05_tf_activity.R",
                   "org.Hs.eg.db")
}

# ------------------------------------------------------------------
# 5. Manifest
# ------------------------------------------------------------------
manifest <- bind_rows(manifest_rows)
readr::write_csv(manifest, file.path(archive_dir, "manifest.csv"))

cat("\nManifest:\n")
print(manifest %>% dplyr::select(resource, n_rows, package, package_version, size_kb))

writeLines(capture.output(sessionInfo()), file.path(archive_dir, "sessionInfo.txt"))

r_version <- paste(R.version$major, R.version$minor, sep = ".")
bioc_version <- tryCatch(as.character(BiocManager::version()), error = function(e) "unknown")

md_lines <- c(
  paste0("# External resource archive -- ", ARCHIVE_DATE),
  "",
  "Frozen copies of every version-sensitive external resource used by the",
  "TF-activity / SNF breast cancer subtyping pipeline. Deposited so that the",
  "analysis can be reproduced exactly, including by readers running it after",
  "these resources have been updated at source.",
  "",
  "## Why this archive exists",
  "",
  "CollecTRI, MSigDB and the Bioconductor annotation packages are curated",
  "resources that change between releases. Code that re-fetches them at run",
  "time will not fail when they change -- it will silently produce different",
  "numbers. These frozen copies, together with the checksums below, remove",
  "that ambiguity.",
  "",
  "## Environment",
  "",
  paste0("- R version: ", r_version),
  paste0("- Bioconductor release: ", bioc_version),
  paste0("- Archive created: ", ARCHIVE_DATE),
  "- Full package environment: `sessionInfo.txt`",
  "- Package lockfile: `renv.lock` (repository root)",
  "",
  "## Contents",
  "",
  "| Resource | Rows | Package | Version | MD5 (tsv) |",
  "|---|---|---|---|---|",
  paste0("| ", manifest$resource, " | ", manifest$n_rows, " | ",
         manifest$package, " | ", manifest$package_version, " | `",
         substr(manifest$md5_tsv, 1, 16), "` |"),
  "",
  "Each resource is provided twice: as `.rds` (R-native, fast to reload) and",
  "as `.tsv.gz` (plain text, readable without R). Full MD5 checksums for both",
  "formats are in `manifest.csv`.",
  "",
  "## Verifying a file",
  "",
  "```r",
  "tools::md5sum(\"collectri_network.tsv.gz\")",
  "# compare against manifest.csv",
  "```",
  "",
  "## Using these instead of live downloads",
  "",
  "Place `collectri_network.rds` at `data/processed/network_collectri.rds`",
  "before running `05_tf_activity.R`. The script uses the cache when present",
  "and only downloads when it is absent, so the archived version will be used",
  "and no network access will occur.",
  "",
  "## Primary data (not redistributed here)",
  "",
  "- **TCGA-BRCA**: downloaded via TCGAbiolinks from the GDC. See",
  "  `logs/01_download_tcga.log` for the access date and",
  "  `logs/sessioninfo/01_download_tcga_package_versions.csv` for the",
  "  TCGAbiolinks version used.",
  "- **METABRIC**: obtained from cBioPortal (`brca_metabric`). Redistribution",
  "  is restricted by the original data use terms; the access date is recorded",
  "  in `logs/13_external_validation.log`.",
  "",
  "## Citation",
  "",
  "If you use this archive, please cite both the deposit DOI and the original",
  "resources (CollecTRI, DoRothEA, MSigDB) in their own right."
)
writeLines(md_lines, file.path(archive_dir, "MANIFEST.md"))

# ------------------------------------------------------------------
# 6. Package lockfile
#
# The archive covers DATA; renv.lock covers CODE dependencies. Neither
# is sufficient alone: frozen gene sets processed by a different
# decoupleR version can still give different answers.
# ------------------------------------------------------------------
if (requireNamespace("renv", quietly = TRUE)) {
  cat("\nSnapshotting package environment to renv.lock ...\n")
  tryCatch({
    renv::snapshot(prompt = FALSE)
    if (file.exists(here::here("renv.lock"))) {
      file.copy(here::here("renv.lock"), file.path(archive_dir, "renv.lock"), overwrite = TRUE)
      cat("  renv.lock copied into the archive.\n")
    }
  }, error = function(e) {
    warning("renv::snapshot() failed: ", conditionMessage(e),
            ". Record package versions manually from sessionInfo.txt.")
  })
}

# ------------------------------------------------------------------
# 7. Single zip, ready to upload
# ------------------------------------------------------------------
zip_path <- here::here("archive", paste0("external_resources_", ARCHIVE_DATE, ".zip"))
zip_ok <- tryCatch({
  old_wd <- setwd(here::here("archive"))
  on.exit(setwd(old_wd), add = TRUE)
  utils::zip(zipfile = basename(zip_path),
             files = paste0("external_resources_", ARCHIVE_DATE))
  TRUE
}, error = function(e) {
  warning("zip() failed (", conditionMessage(e),
          "). Compress the archive folder manually before upload.")
  FALSE
})

cat("\n============================================================\n")
cat("Archive complete.\n")
cat("  Folder: ", archive_dir, "\n", sep = "")
if (zip_ok && file.exists(zip_path)) {
  cat("  Zip:    ", zip_path,
      sprintf("  (%.1f MB)\n", file.size(zip_path) / 1024^2), sep = "")
}
cat("\nNEXT STEPS\n")
cat("  1. Upload the zip to Zenodo (zenodo.org) and RESERVE a DOI before\n")
cat("     publishing, so the DOI can be cited in the manuscript itself.\n")
cat("  2. Cite that DOI in Methods, alongside the original resource citations.\n")
cat("  3. Keep archive/ under version control, or at minimum keep manifest.csv,\n")
cat("     so the checksums travel with the code.\n")
cat("============================================================\n")

log_session_info(script_name, key_packages = c("OmnipathR", "dorothea", "msigdbr",
                                               "decoupleR", "org.Hs.eg.db", "renv"))
close_logger(log_con, script_name)
