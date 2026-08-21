# ============================================================
# 05b_fetch_collectri.R
# Diagnose why decoupleR::get_collectri() failed, and obtain
# CollecTRI by a documented fallback route if it cannot be fixed.
# ============================================================
# CONTEXT
# ------------------------------------------------------------
# 05_tf_activity.R reported:
#
#   CollecTRI could not be retrieved (argument is of length zero).
#
# That message means the call reached OmnipathR but something returned
# empty and a downstream index failed -- it is NOT a "package missing"
# error. The usual causes, in the order worth checking:
#
#   1. OmnipathR or decoupleR too old for each other's API
#   2. Network blocked (university proxy / firewall / TLS interception)
#   3. A corrupted OmnipathR download cache
#
# Run this script. It works through those in order and, only if none
# can be fixed, falls back to the OmniPath web service directly.
#
# WHY THIS MATTERS ENOUGH TO WARRANT ITS OWN SCRIPT
# Without CollecTRI, the grid is DoRothEA-only. Given that the observed
# cross-view concordance between the DoRothEA views is 0.995-1.000 --
# i.e. they are effectively the same view measured twice -- a
# DoRothEA-only grid does not defend against the objection it was built
# to defend against. CollecTRI is currently the single most valuable
# missing piece of the analysis.
#
# OUTPUT: data/processed/network_collectri.rds
# 05_tf_activity.R uses that cache when present and only downloads when
# it is absent, so simply creating it is enough -- no code change.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})
source(here::here("R", "utils.R"))

CACHE_PATH <- here::here("data/processed/network_collectri.rds")

cat("=========================================================\n")
cat("CollecTRI retrieval -- diagnosis and fallback\n")
cat("=========================================================\n\n")

# ------------------------------------------------------------------
# STEP 1. Environment
# ------------------------------------------------------------------
cat("STEP 1: package versions\n")
pkg_version <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(p))
}
versions <- tibble(
  package = c("OmnipathR", "decoupleR", "curl", "httr"),
  version = vapply(c("OmnipathR", "decoupleR", "curl", "httr"), pkg_version, character(1))
)
print(versions)

if (is.na(versions$version[versions$package == "OmnipathR"])) {
  cat("\n  -> OmnipathR is NOT INSTALLED. Install it and re-run:\n")
  cat("     BiocManager::install('OmnipathR')\n\n")
} else {
  # decoupleR's get_collectri() has tracked OmnipathR's API across
  # several releases. A mismatched pair is the most common cause of
  # exactly this error.
  cat("\n  -> If OmnipathR < 3.12 or decoupleR < 2.8, update BOTH together:\n")
  cat("     BiocManager::install(c('OmnipathR','decoupleR'), update = TRUE, ask = FALSE)\n\n")
}

# ------------------------------------------------------------------
# STEP 2. Connectivity
# ------------------------------------------------------------------
cat("STEP 2: can this machine reach omnipathdb.org?\n")
reachable <- FALSE
if (requireNamespace("curl", quietly = TRUE)) {
  h <- curl::new_handle(timeout = 30L)
  probe <- tryCatch(curl::curl_fetch_memory("https://omnipathdb.org/queries/interactions", handle = h),
                    error = function(e) NULL)
  if (!is.null(probe) && probe$status_code == 200) {
    reachable <- TRUE
    cat("  -> reachable (HTTP 200)\n\n")
  } else {
    cat("  -> NOT reachable.",
        if (is.null(probe)) "Connection failed entirely." else paste0("HTTP ", probe$status_code), "\n")
    cat("     On a university/corporate network this is usually a proxy. Try:\n")
    cat("       Sys.setenv(https_proxy = 'http://your.proxy:port')\n")
    cat("     or run this script from a home/unfiltered connection.\n\n")
  }
} else {
  cat("  -> curl not available; skipping connectivity probe.\n\n")
}

# ------------------------------------------------------------------
# STEP 3. Clear a possibly corrupt cache, then retry the proper route
#
# The OmnipathR route is STRONGLY PREFERRED over the fallback in step 4,
# because decoupleR derives the `mor` column from CollecTRI's sign
# annotations by its own documented rules. Reconstructing mor by hand
# risks subtle disagreement with what every other decoupleR user has.
# ------------------------------------------------------------------
cat("STEP 3: retry via OmnipathR (preferred route)\n")
collectri <- NULL

if (requireNamespace("OmnipathR", quietly = TRUE)) {

  # 3a. direct OmnipathR call -- isolates whether the failure is in
  #     OmnipathR or in decoupleR's wrapper around it
  cat("  3a. OmnipathR::collectri() ... ")
  ct_raw <- tryCatch(
    OmnipathR::collectri(organism = 9606L, genesymbols = TRUE, loops = TRUE),
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(ct_raw)) cat("ok --", nrow(ct_raw), "rows\n")

  # 3b. wipe cache and retry, if 3a failed
  if (is.null(ct_raw)) {
    cat("  3b. wiping OmnipathR cache and retrying ... ")
    tryCatch({
      if ("omnipath_cache_wipe" %in% getNamespaceExports("OmnipathR")) {
        OmnipathR::omnipath_cache_wipe()
      }
      ct_raw <- OmnipathR::collectri(organism = 9606L, genesymbols = TRUE, loops = TRUE)
      cat("ok --", nrow(ct_raw), "rows\n")
    }, error = function(e) cat("FAILED:", conditionMessage(e), "\n"))
  }

  # 3c. decoupleR's own wrapper -- what script 05 actually calls
  cat("  3c. decoupleR::get_collectri() ... ")
  ct_dc <- tryCatch(
    decoupleR::get_collectri(organism = "human", split_complexes = FALSE),
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(ct_dc)) {
    cat("ok --", nrow(ct_dc), "rows\n")
    collectri <- ct_dc
  }
}

# ------------------------------------------------------------------
# STEP 4. FALLBACK: OmniPath web service directly
#
# Use ONLY if step 3 cannot be made to work. The `mor` column is
# reconstructed here from is_stimulation / is_inhibition. That follows
# CollecTRI's own sign convention, but it is a REIMPLEMENTATION rather
# than decoupleR's code path, so the two could differ on edge cases.
# If this route is used, say so explicitly in the methods section --
# do not present it as decoupleR::get_collectri() output.
# ------------------------------------------------------------------
if (is.null(collectri) && reachable) {
  cat("\nSTEP 4: FALLBACK -- OmniPath web service directly\n")
  url <- paste0("https://omnipathdb.org/interactions",
                "?datasets=collectri",
                "&genesymbols=yes",
                "&organisms=9606",
                "&fields=sources,references")
  cat("  URL:", url, "\n  downloading ... ")

  raw <- tryCatch(readr::read_tsv(url, show_col_types = FALSE),
                  error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL })

  if (!is.null(raw) && nrow(raw) > 0) {
    cat("ok --", nrow(raw), "rows\n")

    src_col <- intersect(c("source_genesymbol", "source"), colnames(raw))[1]
    tgt_col <- intersect(c("target_genesymbol", "target"), colnames(raw))[1]

    stim <- as.integer(raw[["is_stimulation"]])
    inhi <- as.integer(raw[["is_inhibition"]])

    n_ambiguous <- sum((stim == 1 & inhi == 1) | (stim == 0 & inhi == 0), na.rm = TRUE)
    cat("  edges with ambiguous or absent sign:", n_ambiguous,
        sprintf(" (%.1f%%) -- assigned mor = +1\n", 100 * n_ambiguous / nrow(raw)))

    collectri <- tibble(
      source = raw[[src_col]],
      target = raw[[tgt_col]],
      mor = dplyr::case_when(
        stim == 1 & inhi == 0 ~  1,
        stim == 0 & inhi == 1 ~ -1,
        TRUE                  ~  1     # CollecTRI's default for unsigned edges
      )
    ) %>%
      filter(!is.na(source), !is.na(target), source != "", target != "") %>%
      distinct(source, target, .keep_all = TRUE)

    attr(collectri, "provenance") <- paste0(
      "OmniPath web service fallback (", url, "), retrieved ", Sys.Date(),
      ". mor reconstructed from is_stimulation/is_inhibition, NOT via ",
      "decoupleR::get_collectri(). Report this in the methods section."
    )
    cat("  reconstructed network:", nrow(collectri), "edges,",
        dplyr::n_distinct(collectri$source), "TFs\n")
  }
}

# ------------------------------------------------------------------
# STEP 5. Save
# ------------------------------------------------------------------
cat("\n=========================================================\n")
if (is.null(collectri) || nrow(collectri) == 0) {
  cat("RESULT: CollecTRI could NOT be retrieved by any route.\n\n")
  cat("Do not proceed to the full grid without it. Options, best first:\n")
  cat("  1. Update OmnipathR and decoupleR together, then re-run this script.\n")
  cat("  2. Run this script from an unfiltered network connection.\n")
  cat("  3. Ask a colleague on a working setup to run:\n")
  cat("       saveRDS(decoupleR::get_collectri('human', FALSE), 'collectri.rds')\n")
  cat("     and place the file at data/processed/network_collectri.rds\n\n")
  cat("If none is possible, state in the manuscript that the grid varied\n")
  cat("INFERENCE METHOD and REGULON CONFIDENCE TIER but not REGULON SOURCE,\n")
  cat("and list that as an explicit limitation. Do not leave it unstated.\n")
} else {
  saveRDS(collectri, CACHE_PATH)
  cat("RESULT: CollecTRI cached to\n  ", CACHE_PATH, "\n\n", sep = "")
  cat("  edges:", nrow(collectri), "\n")
  cat("  TFs:  ", dplyr::n_distinct(collectri$source), "\n")
  cat("  mor distribution:\n")
  print(table(collectri$mor))
  prov <- attr(collectri, "provenance")
  if (!is.null(prov)) cat("\n  PROVENANCE:", prov, "\n")
  cat("\nNext: re-run 05_tf_activity.R. It will pick up this cache\n")
  cat("automatically and add the CollecTRI views to the grid.\n")
}
cat("=========================================================\n")
