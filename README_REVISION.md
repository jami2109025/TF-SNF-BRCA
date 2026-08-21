# Revision: from an inconclusive null to a defensible finding

## The problem these scripts solve

Your pipeline's headline result was:

| k | ARI SNF-fused | ARI expression-only | Added value | 95% CI |
|---|---|---|---|---|
| 2 | 0.411 | 0.402 | +0.009 | −0.016 to 0.058 |
| 3 | 0.326 | 0.305 | +0.021 | −0.117 to 0.055 |

Both intervals contain zero, so as reported this supports no claim at all. A CI containing zero is equally compatible with "TF activity adds nothing" (a finding) and "this study could not have detected it" (a non-finding), and nothing in the original design distinguishes them.

Three further problems compound it:

- **The endpoint was rigged.** ARI-vs-PAM50 was the only outcome measure, and PAM50 labels are themselves derived from gene expression. The expression-only baseline is near-optimal for that target by construction. One arm was playing at home.
- **One configuration was tested.** One regulon (DoRothEA A–C), one method (VIPER). "You configured it badly" is an unanswerable objection against a single-cell design.
- **Two large confounds were measured and dropped.** Purity vs PC1 ρ = 0.50 / PC2 ρ = −0.57 (script 04, never followed up), and the METABRIC survival signal collapsing from p = 1.1e-7 to HR 1.09, p = 0.37 once ER was adjusted for (script 13, reported and not discussed).

## Files

```
R/utils_benchmark.R                    new — equivalence testing, cached distances, endpoints
scripts/05_tf_activity.R               REPLACES yours — regulon × method grid
scripts/07_benchmark_models.R          REPLACES yours — grid-wide added value + TOST
scripts/14_prognostic_added_value.R    new — non-circular outcome endpoint (METABRIC)
scripts/15_purity_sensitivity.R        new — the purity confound
scripts/16_incremental_prognostic_value.R  new — value over ER/PAM50; time-varying hazard
```

`R/utils.R` is **not modified**. Its audit trail stays valid; new helpers live in a separate file. Source both:

```r
source(here::here("R", "utils.R"))
source(here::here("R", "utils_benchmark.R"))
```

## The five changes that matter

**1. Equivalence testing (TOST) against a pre-specified SESOI.** This is what makes a negative result publishable. Instead of failing to reject "effect = 0", you test and reject "effect ≥ 0.05 ARI". Rejecting that is positive evidence of absence. Every added-value row now returns one of four verdicts — SUPERIOR, EQUIVALENT, trivially-small, or INCONCLUSIVE — plus the minimum detectable difference, so the power question is answered before a reviewer raises it. SESOI justification is written into `utils_benchmark.R`; do not leave it unargued in the manuscript.

**2. A grid, not a configuration.** DoRothEA A–C, DoRothEA A–B, CollecTRI, plus a degree-preserving **edge-permuted** control for each, crossed with VIPER / ULM / MLM / norm_wsum / consensus. The edge-permuted control is the sharper one: it keeps each TF's out-degree but randomises which genes it regulates, separating "curated regulatory biology helps" from "averaging a few hundred arbitrary gene sets helps". BH correction is applied across views within each (k, endpoint) family — reporting the best of ~15 views unadjusted would be a forking-paths error.

**3. A non-circular endpoint.** Script 14 measures added value as prognostic discrimination of overall survival in METABRIC. Survival is not downstream of the expression matrix, METABRIC took no part in deriving the clusters, and it has ~1,144 events against TCGA's ~150. Neither arm has a structural advantage. **This should be your primary added-value analysis**; the PAM50 comparison becomes a secondary, explicitly circular one.

**4. Cached distances.** d(i,j) is a property of the pair, so it does not change under resampling. The O(n²·p) distance step is computed once and each replicate indexes `D[idx, idx]`, leaving only the cheap affinity/kNN and SNF steps. `verify_distance_cache()` asserts equality against `SNFtool::dist2()` at runtime and **fails hard** rather than producing fast wrong answers, and script 07 additionally rebuilds your primary fused network and compares it to the stored `W_fused.rds`. Your single-view bootstrap took 22.8 hours; the whole grid is now feasible.

**5. Paired resampling, plus a degeneracy guard.** All views see the same patients within a replicate, so the difference is a paired contrast. Separately: your original bootstrap resampled with replacement and had no guard against duplicated patients producing zero distances — `affinityMatrix` divides by the mean kNN distance, so a replicate could return NaN affinities and still produce a finite, meaningless ARI. Now guarded and asserted finite.

## Run order

```
05 → 06 → 07 → (08–13 unchanged) → 14 → 15 → 16
```

Script 05 still writes `tf_activity_viper_AC_primary.rds` with the unchanged primary view, so scripts 06 and 12 run untouched and **every previously reported number stays exactly reproducible**. The revision is additive.

Rough runtimes: 05 ≈ 30–90 min · 07 ≈ 3–8 h for the full grid (vs weeks uncached) · 14 ≈ 30 min · 15 ≈ 20 min · 16 ≈ 10 min.

**Trim the grid while developing.** In `05_tf_activity.R` set `METHODS_REQUESTED <- c("viper","ulm")` and `INCLUDE_EDGEPERM <- FALSE`; in `07_benchmark_models.R` set `N_BOOT <- 50`. Run the full grid once the pipeline is verified end to end.

## Prerequisites

- **CollecTRI** needs `OmnipathR` and internet on first run; it is then cached to `data/processed/network_collectri.rds`. Archive that cache with the manuscript — a live web resource that changes between your run and a reviewer's would silently alter every downstream number.
- Script 15 needs `data/external/tumor_purity_clean.csv` (already present, per your script-04 log).
- Scripts 14/16 need the METABRIC files and, for 16, `table_metabric_assignment_k*.csv` from script 13.

## How to read the output

Read the `conclusion` and `equivalent_after_BH` columns in `table_added_value_equivalence_grid.csv` — **not the CI alone**. Only rows marked EQUIVALENT support a positive claim of absence. Rows marked INCONCLUSIVE support no claim in either direction, and `mdd_80pct_power` tells you the effect size the design could actually have detected.

The grid-level verdict in `table_added_value_grid_verdict.csv` is the sentence your supervisor is asking for.

## What each result would mean

- **All views EQUIVALENT on both endpoints** → your strongest outcome. "TF-activity integration adds nothing beyond expression for BRCA subtyping, demonstrated across a systematically varied grid with equivalence testing on a non-circular endpoint." Methods/benchmarking venue.
- **A view SUPERIOR on the METABRIC prognostic endpoint after BH** → a positive result the PAM50 endpoint was incapable of detecting. Lead with it.
- **All INCONCLUSIVE** → underpowered; make no claim. Check `mdd_80pct_power` against the SESOI and say so plainly.
- **Script 16 T4 interaction survives** → your biological finding: a cluster whose prognostic direction reverses at 10 years, independent of ER. If it does not survive, the crossover is ER's well-known late-recurrence pattern and must not be claimed.

## Three things to be careful about

**Do not report script 13's unadjusted METABRIC log-rank p = 1.14e-07 as external validation.** Your own log shows it collapsing to HR 1.09, p = 0.37 after age + ER. Script 16 replaces it with the nested-model ladder; the decisive row is **M3 vs M2** (clusters on top of age + ER + an existing subtype call), not M3 vs M1.

**The SESOI must be pre-specified and justified in writing.** An unargued equivalence bound is the easiest thing to attack in an equivalence paper. The justification is drafted in `utils_benchmark.R` — put it in the methods section, and do not change the value after seeing results.

**Script 14's bootstrap holds the TCGA centroids fixed.** That is the right frame for evaluating a pre-trained classifier on new data, but it does not propagate discovery-cohort uncertainty. It is stated as a limitation in `NOTE_prognostic_added_value_design.txt`; keep it there rather than letting a reviewer find it.

## Statistical claims these scripts do *not* let you make

Worth being explicit, since the whole point is defensibility:

- Not that SNF fusion is superior to naive concatenation — script 07 tests it, and your existing numbers (0.410 vs 0.411 at k=2) suggest it is not.
- Not that the master regulators in script 12 are discoveries. Your own note flags the circularity, and FOXA1/GATA3/PGR/E2F in BRCA subtypes are textbook. Keep that section as a biological sanity check.
- Not that the clusters are novel subtypes. At k=2 they are essentially the ER+/ER− split.

The contribution is the evaluation framework and what it shows — not a new subtype.
