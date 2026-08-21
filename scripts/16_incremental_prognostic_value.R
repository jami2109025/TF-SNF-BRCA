# ============================================================
# 16_incremental_prognostic_value.R   (NEW)
# Do the SNF clusters add prognostic information beyond ER and PAM50,
# and is the time-varying hazard real or an ER artefact?
# ============================================================
# WHY THIS SCRIPT EXISTS
# ------------------------------------------------------------
# Script 13 reported a striking external-validation result:
#
#     METABRIC k=3, unadjusted log-rank   p = 1.14e-07
#
# and then, in the very same log, the adjusted model that undercuts it:
#
#     assigned_SNF_clusterSNF_C3, adjusted for age + ER:
#         HR = 1.09,  p = 0.37
#
# Almost the entire unadjusted signal is explained by ER status. As it
# stands, the "external validation" mostly demonstrates that the
# clusters re-encode ER -- which is not a finding, because ER status is
# measured routinely by immunohistochemistry and costs nothing.
#
# The manuscript needs to state plainly what survives adjustment. And
# one thing does:
#
#     SNF_C3 x period>10y interaction:  HR = 0.324,  p = 6.97e-11
#
# That is a hazard REVERSAL at ten years -- a group at elevated risk
# early and at reduced risk late -- and it persists after ER
# adjustment. It is the most interesting signal in the entire project
# and is currently buried inside a proportional-hazards caveat, framed
# as a nuisance that invalidates the hazard ratios rather than as the
# result. Non-proportional hazards are only a nuisance when the
# time-varying effect is not the point. Here it is the point.
#
# This script asks the two questions that decide whether any of this is
# publishable as biology:
#
#   Q1  Do the clusters add prognostic information BEYOND age, ER and
#       PAM50/claudin subtype? Nested Cox models, LRT and bootstrapped
#       delta-C-index. If the answer is no, the honest framing is that
#       the clusters recapitulate known subtypes and the contribution is
#       purely methodological.
#
#   Q2  Is the ten-year hazard crossover a property of the CLUSTERS, or
#       is it the well-known ER+ late-recurrence pattern showing
#       through? Tested by fitting the same time-split structure to ER
#       alone, to the clusters alone, and to both together. If ER's own
#       time-split term absorbs the interaction, the crossover is not a
#       cluster finding and must not be claimed as one.
#
# Q2 is the question that decides whether this is a paper or a
# re-description. It is designed here so that it CAN come out negative.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(survival)
  library(survminer)
})
source(here::here("R", "utils.R"))
source(here::here("R", "utils_benchmark.R"))

script_name <- "16_incremental_prognostic_value"
log_con <- init_logger(script_name)
ensure_dirs()
set_pipeline_seed()

N_BOOT <- 2000
CUT_YEARS <- 10          # pre-specified split point, matching script 13
K_VALUES <- c(2, 3)

# ------------------------------------------------------------------
# Inputs: reuse script 13's METABRIC assignments verbatim, so this
# analysis and the published external validation cannot diverge.
# ------------------------------------------------------------------
metabric_clin_path <- here::here("data/external/metabric_clinical.txt")
if (!file.exists(metabric_clin_path)) stop("METABRIC clinical file not found.")

assign_files <- c(k2 = "results/tables/table_metabric_assignment_k2.csv",
                  k3 = "results/tables/table_metabric_assignment_k3.csv")
for (f in assign_files) if (!file.exists(here::here(f))) {
  stop("Missing ", f, " -- run 13_external_validation.R first.")
}
assignments <- lapply(assign_files, function(f) read_csv(here::here(f), show_col_types = FALSE))

metabric_clin <- read.delim(metabric_clin_path, check.names = FALSE, comment.char = "#")

id_col        <- intersect(c("PATIENT_ID", "SAMPLE_ID"), colnames(metabric_clin))[1]
os_months_col <- intersect(c("OS_MONTHS", "os_months"), colnames(metabric_clin))[1]
os_status_col <- intersect(c("OS_STATUS", "os_status"), colnames(metabric_clin))[1]
age_col       <- pick_best_covariate(c("AGE_AT_DIAGNOSIS"), metabric_clin)
er_col        <- pick_best_covariate(c("ER_IHC", "ER_STATUS"), metabric_clin)
subtype_col   <- pick_best_covariate(c("CLAUDIN_SUBTYPE", "THREEGENE", "PAM50"), metabric_clin)
if (any(is.na(c(id_col, os_months_col, os_status_col, age_col, er_col, subtype_col)))) {
  stop("A required METABRIC clinical column could not be located (id/OS/age/ER/",
       "subtype). The subtype column in particular is essential: without it the ",
       "key M3-vs-M2 comparison degenerates into the weaker M3-vs-M1 comparison, ",
       "which would overstate what the clusters add. Columns available: ",
       paste(colnames(metabric_clin), collapse = ", "))
}
cat("Columns used -- id:", id_col, "| age:", age_col, "| ER:", er_col,
    "| subtype:", subtype_col, "\n")

clin <- metabric_clin %>%
  transmute(
    sample_id = .data[[id_col]],
    os_years  = suppressWarnings(as.numeric(.data[[os_months_col]])) / 12,
    os_raw    = as.character(.data[[os_status_col]]),
    age       = suppressWarnings(as.numeric(.data[[age_col]])),
    er_raw    = as.character(.data[[er_col]]),
    subtype   = as.character(.data[[subtype_col]])
  ) %>%
  mutate(
    os_status = case_when(
      grepl("^1|DECEASED|DEAD", os_raw, ignore.case = TRUE) ~ 1,
      grepl("^0|LIVING|ALIVE", os_raw, ignore.case = TRUE) ~ 0,
      TRUE ~ NA_real_),
    er = case_when(
      grepl("^posit", er_raw, ignore.case = TRUE) ~ "Positive",
      grepl("^negat", er_raw, ignore.case = TRUE) ~ "Negative",
      TRUE ~ NA_character_),
    subtype = ifelse(is.na(subtype) | subtype %in% c("", "NC"), NA_character_, subtype)
  )

# ------------------------------------------------------------------
# Q1. NESTED MODELS -- incremental value over routine clinical data
#
# Model ladder, each step adding one block:
#   M0  age
#   M1  age + ER                 (routine IHC -- free, universally available)
#   M2  age + ER + subtype       (expression-based subtype already in the literature)
#   M3  age + ER + subtype + SNF cluster   (this thesis's contribution)
#
# The only number that can justify the clustering clinically is
# M3 vs M2: what the SNF clusters add on top of a subtype call that
# already exists. M3 vs M1 is a weaker and more flattering comparison
# and should not be reported alone.
# ------------------------------------------------------------------
fit_ladder <- function(d) {
  m0 <- coxph(Surv(os_years, os_status) ~ age, data = d)
  m1 <- coxph(Surv(os_years, os_status) ~ age + er, data = d)
  m2 <- coxph(Surv(os_years, os_status) ~ age + er + subtype, data = d)
  m3 <- coxph(Surv(os_years, os_status) ~ age + er + subtype + cl, data = d)
  list(M0 = m0, M1 = m1, M2 = m2, M3 = m3)
}

lrt_step <- function(small, big) {
  stat <- 2 * (as.numeric(logLik(big)) - as.numeric(logLik(small)))
  df   <- length(coef(big)) - length(coef(small))
  tibble(lrt_chisq = stat, df = df, p_value = pchisq(stat, df, lower.tail = FALSE))
}

ladder_results <- list()
boot_results   <- list()

for (k in K_VALUES) {
  asg <- assignments[[paste0("k", k)]] %>%
    dplyr::select(sample_id, assigned_cluster, high_confidence)

  d <- clin %>%
    inner_join(asg, by = "sample_id") %>%
    filter(is.finite(os_years), os_years > 0, !is.na(os_status),
           !is.na(age), !is.na(er), !is.na(subtype), !is.na(assigned_cluster)) %>%
    mutate(cl = factor(assigned_cluster), er = factor(er), subtype = factor(subtype))

  cat("\n========== k =", k, "==========\n")
  cat("Complete-case N:", nrow(d), " events:", sum(d$os_status), "\n")
  cat("Subtype levels:", paste(levels(d$subtype), collapse = ", "), "\n")

  m <- fit_ladder(d)

  step_table <- bind_rows(
    bind_cols(tibble(step = "M1 vs M0: + ER"),                        lrt_step(m$M0, m$M1)),
    bind_cols(tibble(step = "M2 vs M1: + subtype"),                   lrt_step(m$M1, m$M2)),
    bind_cols(tibble(step = "M3 vs M2: + SNF cluster (KEY TEST)"),    lrt_step(m$M2, m$M3)),
    bind_cols(tibble(step = "M3 vs M1: + subtype + SNF cluster"),     lrt_step(m$M1, m$M3))
  ) %>% mutate(k = k, .before = 1)

  cindex_table <- tibble(
    k = k,
    model = c("M0 age", "M1 age+ER", "M2 age+ER+subtype", "M3 age+ER+subtype+cluster"),
    cindex = c(summary(m$M0)$concordance[1], summary(m$M1)$concordance[1],
               summary(m$M2)$concordance[1], summary(m$M3)$concordance[1])
  )

  cat("\nNested model LRTs:\n"); print(step_table, width = Inf)
  cat("\nC-index ladder:\n"); print(cindex_table)

  ladder_results[[paste0("k", k)]] <- list(steps = step_table, cindex = cindex_table)

  # Bootstrap the incremental C-index of M3 over M2. The LRT alone is
  # not enough: with ~1900 patients and ~1100 events, a trivially small
  # improvement can be highly significant. The delta-C-index says
  # whether it is also large enough to matter.
  cat("\nBootstrapping delta C-index (M3 - M2), n =", N_BOOT, "...\n")
  rows <- seq_len(nrow(d))
  d_cidx <- numeric(N_BOOT)
  for (b in seq_len(N_BOOT)) {
    set.seed(PIPELINE_SEED + 7000 + b)
    db <- d[sample(rows, length(rows), replace = TRUE), ]
    fit2 <- try(coxph(Surv(os_years, os_status) ~ age + er + subtype, data = db), silent = TRUE)
    fit3 <- try(coxph(Surv(os_years, os_status) ~ age + er + subtype + cl, data = db), silent = TRUE)
    d_cidx[b] <- if (inherits(fit2, "try-error") || inherits(fit3, "try-error")) NA_real_ else
      summary(fit3)$concordance[1] - summary(fit2)$concordance[1]
  }
  eq <- tost_equivalence_bootstrap(d_cidx, sesoi = 0.02,
                                   label = paste0("delta_cindex_M3_over_M2_k", k))
  eq$k <- k
  boot_results[[paste0("k", k)]] <- eq
  cat("\nDelta C-index (SNF clusters over age+ER+subtype):\n"); print(eq, width = Inf)
}

steps_all  <- bind_rows(lapply(ladder_results, `[[`, "steps"))
cindex_all <- bind_rows(lapply(ladder_results, `[[`, "cindex"))
boot_all   <- bind_rows(boot_results)

write_csv(steps_all,  here::here("results/tables/table_metabric_nested_model_lrt.csv"))
write_csv(cindex_all, here::here("results/tables/table_metabric_cindex_ladder.csv"))
write_csv(boot_all,   here::here("results/tables/table_metabric_incremental_cindex_bootstrap.csv"))

# ------------------------------------------------------------------
# Q2. TIME-VARYING EFFECT -- cluster-specific, or ER showing through?
#
# Four models on the same episode-split data:
#   T1  cluster x period                    (is there a crossover at all?)
#   T2  ER x period                         (does ER alone show the same?)
#   T3  cluster x period + ER               (crossover with ER held constant)
#   T4  cluster x period + ER x period      (BOTH allowed to vary in time --
#                                            the decisive test)
#
# T4 is the one that matters. Only if the cluster x period interaction
# survives while ER is ALSO permitted its own time-varying effect can
# the crossover be attributed to the clusters. Fitting only T1 or T3
# and declaring a novel time-varying subtype effect would be the
# classic error here, because ER's own well-documented late-recurrence
# pattern would be forced into the cluster term.
# ------------------------------------------------------------------
timevar_results <- list()
rmst_results <- list()

for (k in K_VALUES) {
  asg <- assignments[[paste0("k", k)]] %>% dplyr::select(sample_id, assigned_cluster)
  d <- clin %>%
    inner_join(asg, by = "sample_id") %>%
    filter(is.finite(os_years), os_years > 0, !is.na(os_status), !is.na(age), !is.na(er)) %>%
    mutate(cl = factor(assigned_cluster), er = factor(er)) %>%
    dplyr::select(-os_raw, -er_raw)   # character columns survSplit does not need

  # NOTE: the id variable is named `subject_id`, not `id`. survSplit()
  # CREATES the id column, and passing a name that already exists in
  # `data` silently overwrites it -- which would corrupt the per-subject
  # clustering of episodes and hence every robust standard error below.
  dsplit <- survSplit(Surv(os_years, os_status) ~ ., data = d,
                      cut = CUT_YEARS, episode = "period_id", id = "subject_id") %>%
    mutate(period = factor(ifelse(period_id == 1, "early_0_10y", "late_gt10y"),
                           levels = c("early_0_10y", "late_gt10y")))

  T1 <- coxph(Surv(tstart, os_years, os_status) ~ cl * period + age, data = dsplit)
  T2 <- coxph(Surv(tstart, os_years, os_status) ~ er * period + age, data = dsplit)
  T3 <- coxph(Surv(tstart, os_years, os_status) ~ cl * period + er + age, data = dsplit)
  T4 <- coxph(Surv(tstart, os_years, os_status) ~ cl * period + er * period + age, data = dsplit)

  tidy_model <- function(fit, nm) {
    broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(model = nm, k = k, .before = 1)
  }
  tv <- bind_rows(
    tidy_model(T1, "T1_cluster_x_period"),
    tidy_model(T2, "T2_ER_x_period"),
    tidy_model(T3, "T3_cluster_x_period_plus_ER"),
    tidy_model(T4, "T4_cluster_x_period_plus_ER_x_period")
  )
  timevar_results[[paste0("k", k)]] <- tv

  cat("\n===== Time-varying models, k =", k, "=====\n")
  cat("\nT4 (decisive: BOTH cluster and ER allowed time-varying effects):\n")
  print(T4 %>% broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
          dplyr::select(term, estimate, conf.low, conf.high, p.value), n = 30)

  key_terms <- tv %>%
    filter(model == "T4_cluster_x_period_plus_ER_x_period",
           grepl("^cl.*:period|^period.*:cl", term))
  survives <- any(key_terms$p.value < 0.05, na.rm = TRUE)
  cat("\nVERDICT k=", k, ": cluster x period interaction ",
      if (survives) "SURVIVES" else "does NOT survive",
      " with ER also allowed a time-varying effect.\n",
      if (survives)
        "  -> the ten-year hazard crossover is NOT reducible to ER status; report it as a cluster finding.\n"
      else
        "  -> the crossover is attributable to ER; do NOT claim it as a novel cluster property.\n",
      sep = "")

  # PH-assumption-free corroboration: RMST within each period. A Cox
  # interaction term is itself a model-based quantity, so a
  # model-free statistic in the same direction materially strengthens
  # the claim (and would expose it if the interaction were an artefact
  # of the PH parameterisation).
  if (requireNamespace("survRM2", quietly = TRUE) && length(unique(d$cl)) >= 2) {
    for (per in c("early", "late")) {
      dd <- if (per == "early") {
        d %>% mutate(t = pmin(os_years, CUT_YEARS),
                     s = ifelse(os_years <= CUT_YEARS, os_status, 0))
      } else {
        d %>% filter(os_years > CUT_YEARS) %>%
          mutate(t = os_years - CUT_YEARS, s = os_status)
      }
      lv <- levels(droplevels(dd$cl))
      if (length(lv) < 2 || nrow(dd) < 50) next
      for (i in 1:(length(lv) - 1)) for (j in (i + 1):length(lv)) {
        sub <- dd %>% filter(cl %in% c(lv[i], lv[j]))
        arm <- as.numeric(sub$cl == lv[j])
        # tau must lie inside the observed follow-up of BOTH arms, or
        # RMST is extrapolating past the data in one of them. Guard
        # against an arm that is empty or has no positive follow-up
        # after the period restriction -- common in the late window.
        if (sum(arm == 0) < 10 || sum(arm == 1) < 10) next
        tau <- min(max(sub$t[arm == 0]), max(sub$t[arm == 1])) * 0.95
        if (!is.finite(tau) || tau <= 0) next
        r <- try(survRM2::rmst2(sub$t, sub$s, arm, tau = tau), silent = TRUE)
        if (inherits(r, "try-error")) next
        rmst_results[[length(rmst_results) + 1]] <- tibble(
          k = k, period = per, comparison = paste(lv[j], "minus", lv[i]),
          tau = tau,
          rmst_diff = r$unadjusted.result[1, 1],
          ci_lower = r$unadjusted.result[1, 2],
          ci_upper = r$unadjusted.result[1, 3],
          p_value = r$unadjusted.result[1, 4]
        )
      }
    }
  }
}

timevar_all <- bind_rows(timevar_results)
write_csv(timevar_all, here::here("results/tables/table_metabric_timevarying_models.csv"))

if (length(rmst_results) > 0) {
  rmst_all <- bind_rows(rmst_results) %>%
    group_by(k, period) %>% mutate(FDR = p.adjust(p_value, "BH")) %>% ungroup()
  cat("\n=== RMST by period (PH-assumption-free corroboration) ===\n")
  print(rmst_all, width = Inf)
  write_csv(rmst_all, here::here("results/tables/table_metabric_rmst_by_period.csv"))
}

# ------------------------------------------------------------------
# Figure: time-varying hazard ratios, cluster vs ER side by side
# ------------------------------------------------------------------
plot_dat <- timevar_all %>%
  filter(model == "T4_cluster_x_period_plus_ER_x_period",
         grepl("^cl|^er", term), !grepl("^age", term)) %>%
  mutate(kind = ifelse(grepl("^cl", term), "SNF cluster", "ER status"),
         interaction = grepl(":", term))

if (nrow(plot_dat) > 0) {
  p <- ggplot(plot_dat, aes(x = term, y = estimate, colour = kind)) +
    geom_hline(yintercept = 1, colour = "grey40") +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    coord_flip() + scale_y_log10() +
    facet_wrap(~ paste0("k = ", k), scales = "free_y") +
    labs(x = NULL, y = "Hazard ratio (log scale)",
         title = "Time-varying effects in METABRIC: SNF cluster vs ER status",
         subtitle = paste0("Model T4: both cluster and ER allowed period-specific effects (split at ",
                           CUT_YEARS, " years). Interaction terms are the crossover.")) +
    theme_bw(base_size = 10)
  ggsave(here::here("results/figures/figure_metabric_timevarying_cluster_vs_ER.png"),
         p, width = 11, height = 6, dpi = 150)
}

writeLines(c(
  "INCREMENTAL PROGNOSTIC VALUE AND THE TIME-VARYING EFFECT.",
  "======================================================================",
  "",
  "CONTEXT. Script 13 reported METABRIC k=3 log-rank p = 1.14e-07, but",
  "the same log shows SNF_C3 adjusted for age + ER at HR = 1.09,",
  "p = 0.37. Most of the unadjusted signal is ER status. Reporting the",
  "unadjusted p-value as external validation without that adjustment",
  "would be misleading, because ER is measured routinely by IHC and",
  "costs nothing.",
  "",
  "Q1 -- WHAT THE CLUSTERS ADD [table_metabric_nested_model_lrt.csv,",
  "table_metabric_cindex_ladder.csv,",
  "table_metabric_incremental_cindex_bootstrap.csv].",
  "The decisive comparison is M3 vs M2: SNF clusters added on top of",
  "age + ER + an EXISTING expression-based subtype call. M3 vs M1",
  "(without subtype) is a weaker and more flattering comparison and",
  "must not be reported on its own. With ~1900 patients and ~1100",
  "events, a trivially small improvement can be highly significant, so",
  "the bootstrapped delta C-index and its equivalence test -- not the",
  "LRT p-value -- determine whether the gain matters.",
  "",
  "Q2 -- IS THE CROSSOVER REAL [table_metabric_timevarying_models.csv].",
  "Script 13 found a strong SNF_C3 x period>10y interaction",
  "(HR = 0.324, p = 7e-11): a hazard REVERSAL at ten years. It is the",
  "most interesting signal in the project and is currently buried in a",
  "proportional-hazards caveat. Non-proportional hazards are a nuisance",
  "only when the time-varying effect is not the point; here it is.",
  "",
  "But ER+ disease is already known to show late recurrence, so the",
  "crossover must be tested against ER rather than merely adjusted for",
  "it. Model T4 lets BOTH the clusters and ER have period-specific",
  "effects. Only a cluster x period interaction that survives in T4 may",
  "be claimed as a cluster property. Fitting only T1 or T3 would force",
  "ER's own late-recurrence pattern into the cluster term -- the",
  "classic error in this exact situation.",
  "",
  "RMST by period is reported as model-free corroboration, because a",
  "Cox interaction is itself a model-based quantity.",
  "",
  "HOW TO WRITE THIS UP.",
  "If T4 survives: this is the manuscript's biological result -- a",
  "cluster whose prognostic direction reverses at ten years,",
  "independent of ER. If it does not: state that the external",
  "validation reflects ER status, and rest the contribution on the",
  "benchmarking framework instead. Both are defensible. Reporting the",
  "unadjusted log-rank p-value as if it were the finding is not.",
  "",
  "LIMITATION. The ten-year split point is pre-specified to match",
  "script 13, not data-driven. A data-chosen cut point would require",
  "correction for the selection, and is not attempted here."
), here::here("results/tables/NOTE_incremental_prognostic_value.txt"))

log_session_info(script_name, key_packages = c("survival", "survRM2", "broom"))
cat("\n✓ 16_incremental_prognostic_value.R complete.\n")

close_logger(log_con, script_name)
