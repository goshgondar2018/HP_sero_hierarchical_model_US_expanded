# NEWS

## Unreleased

- **Switched to unweighted NHANES data, via a new `data/re_extracting_nhanes_data.R`**
  that re-aggregates raw NHANES III (1988-1994) and NHANES 1999-2000 records
  by age/sex/race using unweighted counts (`total`/`HP`), replacing the
  fractional, survey-weighted pseudo-counts previously baked into
  `df_hp_usa_w_final.Rdata` (identified during the earlier review as
  producing pooled race-specific prevalence that didn't match the published
  McMahon et al. Table 1 NHANES estimates, especially for NHW). Output is
  `data/df_hp_usa_w_final_unweighted_nhanes.Rdata`.
  All five model/figures/LOO scripts were updated to `load()` this file in
  place of the old one, but every one of them had the same bug: they loaded
  `data/df_hp_usa_final_unweighted_nhanes.Rdata` (object name
  `df_hp_usa_final`, no "w") while every line downstream still referenced
  `df_hp_usa_w_final` (with "w") -- a plain object-name mismatch. In a fresh
  session this errors out immediately ("object 'df_hp_usa_w_final' not
  found"); in an interactive session with a stale `df_hp_usa_w_final` left
  over from an earlier `load()` of the old file, it would have **silently
  kept using the old weighted data with no error at all**. Fixed by pointing
  all five scripts at the correctly-named
  `df_hp_usa_w_final_unweighted_nhanes.Rdata` (content-identical, verified
  via `all.equal()`), and by fixing `re_extracting_nhanes_data.R` itself to
  save under the `df_hp_usa_w_final` object name directly, so a future
  re-run (e.g. with updated raw NHANES extracts) doesn't reintroduce the same
  mismatch. **Requires refitting** every model, since this changes the
  likelihood input data, not just post-processing.
- **`figures_cluster_sensitivity_analysis_sep15.R` read fit objects from the
  wrong directory.** All three `readRDS()` calls (`base_continuous_bc`,
  `RW_continuous_bc`, `RW_expanded_continuous_bc`) pointed at
  `./results_cluster/`, which only ever receives the 6 *binary*-cohort model
  fits from `all_models_cluster.R`; the continuous-cohort fits are saved to
  `./results_cluster_sensitivity_analysis/` by
  `all_models_sensitivity_analysis_cluster.R`. Fixed to the correct directory
  for all three. Found while adding the effective-sample-size calculation
  below, which reads one of these same fit objects.
- **Added overdispersion-adjusted effective sample size**, alongside the
  existing raw `total_pop`. For `y ~ BetaBinomial(n, p*theta, (1-p)*theta)`,
  `Var(y)` is inflated relative to a plain `Binomial(n, p)` by a design
  effect `DEFF = (n + theta) / (1 + theta)`; the effective sample size is the
  plain-binomial `n` giving the same precision:
  `n_eff = n*(1 + theta) / (n + theta)`, computed per row and per posterior
  draw of `theta`, then summed across rows and reported with a 95% CrI.
  `figures_cluster_sensitivity_analysis_sep15.R` computes this for the
  selected continuous-cohort model only (`RW_continuous_bc`), written to
  `results/effective_sample_size_RWcohort.csv`.
  `figures_cluster_sep4.R` computes this for every model in the binary-cohort
  registry (not just the selected `study_re`), since `total_pop` is the same
  across models but `theta` -- and so `eff_n` -- varies by model; written to
  `results/effective_sample_size_ALL.csv` with one row per model. Requires
  the corresponding model(s) to already be fit (reads the saved `.rds`); no
  refit needed for this calculation itself.

## Unreleased (branch `fix/foi-and-denominators`)

A second review pass, run independently against the same co-author workbook
and the Supplementary Appendix, layered on top of the
`bugfix/data-and-model-review` fixes (recorded in the section below). The
items here are the ones that review did not cover.

**Results will change.** The force-of-infection and denominator corrections
both alter reported numbers. All models must be refitted and all results
regenerated.

### Force of infection was understated by up to 27%

`figures_cluster_sep4.R` and `figures_cluster_sensitivity_analysis_sep15.R`
both computed the force of infection with the denominator

```r
foi_den <- 1 - p_mat * one_minus_exp + 1e-10
```

`p_mat` already **is** F(a), the cumulative probability of infection (the
rescaled `p_raw_mat`). Multiplying it again by `one_minus_exp` applied the
cumulative-infection term a second time, so the denominator was always larger
than the `1 - F(a)` that Appendix Eq. 3 requires and every force-of-infection
estimate was too small.

The bias is age-dependent, and largest where the force of infection peaks
(alpha = 0.8, gamma = 0.05, kappa = 1):

| Age | Correct lambda(a) | As computed | Error |
|----:|------------------:|------------:|------:|
|   5 |           0.03785 |     0.03243 | -14.3% |
|  20 |           0.02977 |     0.02163 | **-27.3%** |
|  40 |           0.01756 |     0.01346 | -23.4% |
|  90 |           0.00213 |     0.00204 |  -4.3% |

Because the bias varies with age, it distorted both the peak force of
infection and the age at which it peaks, and did so differently for each
race/ethnicity group. Recomputing from a short test fit of
`kappa_partial_study_re` gives an indication of the scale (parameter values
from a real run will differ, but the direction and rough magnitude hold):

| Group | Cohort | Peak FOI (old) | Peak FOI (corrected) | Change | Age at peak (old -> corrected) |
|---|---|---:|---:|---:|---|
| Hispanic | Pre-1960 | 0.1077 | 0.1382 | +28.3% | 10.3 -> 14.6 |
| Hispanic | 1960+    | 0.0777 | 0.0922 | +18.7% | 9.5 -> 15.1 |
| AI/AN    | Pre-1960 | 0.1115 | 0.1298 | +16.4% | 2.2 -> 3.1 |
| AI/AN    | 1960+    | 0.1685 | 0.1966 | +16.6% | 1.1 -> 1.3 |
| API      | Pre-1960 | 0.0361 | 0.0410 | +13.6% | 14.0 -> 15.8 |
| API      | 1960+    | 0.0336 | 0.0364 |  +8.3% | 17.0 -> 19.4 |
| NHB      | Pre-1960 | 0.0814 | 0.0878 |  +7.9% | 0.7 -> 0.9 |
| NHB      | 1960+    | 0.0627 | 0.0663 |  +5.7% | 0.7 -> 0.8 |
| NHW      | Pre-1960 | 0.0174 | 0.0190 |  +9.2% | 19.5 -> 20.8 |
| NHW      | 1960+    | 0.0193 | 0.0202 |  +4.3% | 7.2 -> 7.9 |

Peak FOI rises by 4-28% and the age at peak moves **later**, by up to 5.6
years. Because the shift differs by group, between-group FOI comparisons
change as well, not just the levels.

The denominator now uses `p_raw_mat`, matching the prevalence curve reported
from the same loop. Age 0 is excluded from the FOI in both scripts (it was
already excluded in the continuous-cohort script but commented out in the
binary one): `(gamma*a)^(kappa-1)` is infinite at age 0 when kappa < 1.

### Three observations used the wrong binomial denominator

`subgroup_n` (the number enrolled) was passed to JAGS as the binomial
denominator, but `total` (== `hp_denom`, the number actually assayed) is what
the dataset's own `prev` column uses. They differ for three rows:

| Study | Group | Numerator | Old denominator | Correct | Old prevalence | Correct | Error |
|---|---|---:|---:|---:|---:|---:|---:|
| MeiSocS2016 | Hispanic | 1,149 | 1,562 | 1,226 | 73.6% | 93.7% | **20.2 pp** |
| RubGenE2013 | AI/AN | 594 | 801 | 766 | 74.2% | 77.5% | 3.4 pp |
| MieHeli2018 | AI/AN | 422 | 571 | 561 | 73.9% | 75.2% | 1.3 pp |

MeiSocS2016 is a 1,226-person Hispanic observation at age 70 (pre-1960
cohort), so the old denominator biased Hispanic maximum prevalence downward.
Changed in both model scripts, the LOO script and the binary figures script.

### Index of disparity

The Pearcy & Keppel index of disparity was computed as
`(sum(Mean - Mean[5]) / 5) * 100 / Mean[5]`. Three problems, all fixed:

- **`abs()` was missing**, so groups above and below the referent cancelled.
- **The sum was divided by 5** (all groups) rather than by the 4 comparison
  groups; the referent contributes a zero deviation, shrinking the index by a
  further factor of 4/5.
- **The interval was built from marginal quantiles** (`LB - LB[5]`,
  `UB - UB[5]`, then divided by `UB[5]`/`LB[5]`), which is not a valid
  credible interval for a nonlinear function of correlated parameters.

The index is now computed per posterior draw from the existing
`max_prev_draws_pre` / `max_prev_draws_post` matrices, then summarised.

### Change in maximum prevalence is now computed per draw

`max_prev_percent_change_summarized` differenced the already-summarised
`Mean`, `LB` and `UB` columns. Now computed per draw, and reported as both an
absolute change (percentage points) and a relative change (%).

### Analytic sample size, and the provenance of "N = 62,688"

`total_pop = sum(unique(data_model_filtered$subgroup_n))` de-duplicated by
VALUE across the whole dataset, silently dropping every subgroup whose size
coincided with another's -- 22 subgroups. That expression is the source of the
**N = 62,688** reported in the manuscript: it reproduces that figure exactly
from the figure script's dataset, whose correct sum is 64,897. It is not a
correct participant count for that or any other dataset. Now
`sum(data_model_filtered$total)`.

The study and participant counts in the manuscript should be re-derived from
the corrected pipeline before the next submission.

### Two further precision-vs-SD prior mismatches — resolved by author

Following the same pattern as the gamma-intercept prior resolved above, JAGS
`dnorm`/`dt` take a PRECISION, and two more priors were written as though the
second argument were an SD or a variance. Both are named constants, passed to
JAGS as data (`KAPPA_TAU` in `all_models_cluster.R`, `RW_INIT_TAU` in
`all_models_sensitivity_analysis_cluster.R`), with the reasoning documented
inline at each definition. Both now resolved:

- **`KAPPA_TAU`: `0.40` → `1.5`.** As originally coded (0.40 read as
  precision), kappa's 95% prior interval was ~(0.05, 22); the Appendix's
  "kappa ~ Lognormal(0, 0.4)" is itself ambiguous between 0.4 as an SD
  (precision 6.25, kappa ~(0.46, 2.19)) and as a variance (precision 2.5,
  kappa ~(0.29, 3.44)). Author decision: rather than adopt either literal
  reading, target kappa's 95% prior interval directly at ~(0.2, 5)
  (`KAPPA_TAU = 1/(ln(5)/1.96)^2`, rounded to 1.5) — tighter than the
  original but looser than either Appendix reading. Median stays at
  kappa = 1 regardless of `KAPPA_TAU`.
- **`RW_INIT_TAU`: `25` → `0.04`.** As originally coded (25 read as
  precision), the random walk's earliest-cohort state had SD 0.2 — a fairly
  tight prior. The Appendix's "N(0, 25)" reads more naturally as variance =
  25 (SD = 5, precision = 0.04). Author decision: adopt the Appendix's
  diffuse reading, so the earliest birth cohort's random-walk state is left
  sufficiently noninformative.

**Both require refitting** `all_models_cluster.R` and
`all_models_sensitivity_analysis_cluster.R`.

### Convergence diagnostics

- **Chains started from identical values.** The inits set only `.RNG.seed`, so
  all three chains began from JAGS's default initialisation. Gelman-Rubin PSRF
  assumes chains start over-dispersed relative to the posterior; with
  identical starts it is anti-conservative and can indicate convergence that
  has not occurred. The binary-cohort chains now start from dispersed `mub`
  and `log_theta` values. (The continuous models disperse only `log_theta`:
  their birth-cohort slopes multiply `bc_scaled`, which spans ~91 years, so a
  slope start of even 0.5 gives `log(gamma) = -3.9 + 45` and overflows the
  chain. The same scale mismatch affects the Cauchy(0,1) slope PRIOR, which on
  this scale puts most of its mass on implausibly large effects -- centring
  and scaling `bc_scaled`, for example per decade, is worth considering.)

  **Update:** with the original spread (`mub` offsets of ±1.5/0.5/1.0/0.5,
  `log_theta` ±1), running the full 500,000-iteration burn-in in
  `all_models_cluster.R` produced `max_PSRF_all` in the 1.2–2.1 range across
  all 6 models -- the chains weren't merging back together within budget. This
  isn't necessarily a new problem: with identical starts there was no way for
  the diagnostic to detect it in the first place (all chains inherit the same
  bias from the same starting point, so between-chain disagreement -- the only
  thing PSRF measures -- was structurally impossible to observe). Rather than
  lengthen the already-expensive burn-in, the dispersion was halved (author
  decision) -- still genuinely non-identical starts, just a shorter distance
  for the chains to travel back to consensus:
  `mub` offsets ±0.75/0.25/0.5/0.25, `log_theta` ±0.5. **Requires refitting**;
  worth checking whether `max_PSRF_all` comes down to a reasonable range
  before deciding if a longer burn-in is also needed.

  **Continuous-cohort models also don't converge** (`max_PSRF_all`:
  `base_continuous_bc` 1.30, `RW_continuous_bc` 1.18,
  `RW_expanded_continuous_bc` 2.42), despite only `log_theta` being
  deliberately dispersed there. Two distinct causes, not one: (1)
  `base_continuous_bc`'s slope terms multiply `bc_scaled` directly (see
  above) -- a badly-scaled covariate/prior combination that makes for
  inherently slow mixing regardless of initialization. **Author decision:
  left as-is for now** -- rescaling `bc_scaled` (e.g. to decades) would
  likely help, but touches both the model and figures scripts and changes
  what the slope coefficients mean (per-year vs. per-decade), so revisit
  deliberately rather than as a side effect of chasing convergence.
  (2) Trace plots for `RW_continuous_bc` (`par_alpha0`/`par_gamma0` by
  race, `rw_alpha_pop`/`rw_gamma_pop` by cohort) showed the same slow,
  synchronized wandering in every panel regardless of race or cohort --
  chain 2 riding high and chain 3 riding low at the same iterations across
  every one of them. That's the signature of a genuine non-identifiability,
  not just high dimensionality: `par_alpha0[j]`/`mub[1]` and
  `rw_alpha_pop`'s own level enter the likelihood only as a sum
  (`logit(alpha[i]) = par_alpha0[inds[i]] + rw_alpha_pop[cohort_idx[i]] +
  study_re[...]`), so shifting every race intercept up by δ while shifting
  the whole `rw_alpha_pop` trajectory down by δ leaves every prediction
  unchanged -- a flat ridge in the posterior that MCMC wanders along
  slowly. The diffuse `RW_INIT_TAU` (0.04, SD 5) adopted above made this
  worse, not better, by giving `rw_alpha_pop[1]` -- one end of that ridge
  -- more room to wander. **Author decision: revert `RW_INIT_TAU` to 25.00
  (SD 0.2)** as a pragmatic mitigation -- it anchors one end of the ridge
  and constrains the wandering, without removing the underlying
  non-identifiability (a proper fix would need a structural constraint,
  e.g. a sum-to-zero constraint on `rw_alpha_pop`, which is a larger
  reparameterization left for later if this isn't sufficient). This is a
  reasonable stopping point here because every downstream use of these
  parameters (predicted prevalence, FOI, median age of infection) depends
  only on the *sum* `par_alpha0[inds[i]] + rw_alpha_pop[cohort_idx[i]]`,
  which can be well-identified even when the individual components aren't
  -- worth checking `max_PSRF_all` on the derived prevalence curves
  specifically, not just the raw parameters, before concluding more work
  is needed. **Requires refitting.** `RW_expanded_continuous_bc` likely
  has the same underlying issue (same additive structure, plus its own
  race-specific deviations), and remains the least converged of the three;
  not separately addressed here.
- **Convergence was summarised by the mean PSRF**, which hides a single badly
  mixing parameter among 15-40. Both model scripts now also record
  `max_PSRF_all`, `n_PSRF_above_1.1` and a `converged` flag, print the maximum
  alongside the mean, and emit a loud `*** NOT CONVERGED ***` marker and an R
  warning when the maximum exceeds 1.1.

### Further bug fixes

- **`LOO_all_model_comparison.R` could not run**: it uses
  `as.matrix(mcmc_list)` on coda `mcmc.list` objects but never loaded `coda`,
  so `as.matrix` fell through to `as.matrix.default` and returned an unnamed
  list-matrix. Every subsequent `mat_all[, "par[1,1]"]` lookup then failed with
  `no 'dimnames' attribute for array`. Added `library(coda)`.
- **`gtsave()` to PNG aborted the whole figure loop.** It requires `webshot2`,
  `chromote` and a headless Chrome; where those are absent it throws, which
  killed `figures_cluster_sep4.R` after the first model and lost every
  subsequent model's figures for the sake of a diagnostics table. PSRF tables
  are now always written as CSV, with the PNG attempted best-effort and a note
  if it is skipped.
- **Observed prevalence points were unweighted means.** `df_p_curves_obs`
  averaged subgroup prevalences within each race x era x age cell, so a
  15-person subgroup counted as much as a 6,000-person one, and it averaged
  confidence limits (not a valid interval for the pooled cell). Cells are now
  pooled by counts with a Jeffreys binomial interval.
- `summarize_max_prev_RD()` labelled its output `type = "Risk Ratio"` although
  it computes a difference; now "Risk Difference".

### Open questions for the authors (no code change)

- **`TsaAmJE2005` uses combined-modality testing** (ELISA for ages >= 2, stool
  HpSA below 2) and is included, although the manuscript states that
  combined-modality testing was excluded.
- **`MuhJInf2015` has no Study Details record** in the workbook, and
  `ParJInf1997` has a Study Details record but no data rows -- an internal
  inconsistency in the source workbook. `MuhJInf2015` is dropped by the
  current filters, but the mismatch is worth raising with the co-author.
- **NHANES cells can be very small after rounding.** The 250 NHANES rows carry
  rescaled survey-weighted counts; 18 cells have an effective denominator below
  5, and rounding to integers shifts cell prevalence by more than 5 percentage
  points in 13 cells (worst case 54 pp, in an extreme-age cell). A minimum
  effective sample size for these cells is worth considering.
- **The random walk treats unequally spaced cohorts as adjacent.** The 69
  distinct birth-cohort years span 91 calendar years, with 14 steps of 2-4
  years. A unit-step random walk applies the same innovation variance to a
  1-year and a 4-year gap, over-smoothing across the gaps. Scaling each
  innovation precision by the actual gap is the discrete analogue of Brownian
  motion.
- **DIC and ELPD-LOO are not on the same scale.** DIC comes from the JAGS
  `deviance` node, which conditions on the latent `pi[i]`, whereas the PSIS-LOO
  log-likelihood integrates `pi[i]` out analytically. The two rankings are not
  directly comparable; ELPD-LOO is the appropriate basis for model comparison.
- **Manuscript/appendix inconsistencies** noticed while cross-checking: the
  FOI results paragraph cites "(Figure 5)" mid-sentence where the surrounding
  text refers to Figure 3; the group-specific kappa means (1.13, 0.849, 1.73,
  0.894, 1.92) are attributed in the text to the continuous-cohort model but
  listed in Table S1 under the binary-cohort group-specific kappa model with
  study random effects; the Discussion reports the largest maximum-prevalence
  gap widening "from 33 to 58 percentage points" while the Results give
  30.2 -> 60.4 (Hispanic - NHW) and 31.3 -> 56.9 (AI/AN - NHW); equation
  number (23) is used twice; and "(5% CrI: 30.4, 76.7)" is missing a digit.


## Unreleased (branch `bugfix/data-and-model-review`)

A full review of the data and code behind "Widening racial and ethnic
disparities in *Helicobacter pylori* acquisition in the USA, 1965–2014:
a Bayesian meta-analytic catalytic model of age-specific forces of
infection", checked against a newer, co-author-curated data export
(`data/HPintheUS_WorkingFile_V12_20230727.xlsx`). See individual commits
for full detail; summarized here.

### Data fixes

- **`data/df_hp_usa_w_final.Rdata` corrected to match the co-author's
  curated inclusion flags** (`viz_flag_race_eth_seer_comb` in the Excel
  file's "Long HP data" sheet). Three studies had rows that disagreed
  with the flag:
  - **ButGast2019 (Butt et al.)**: all four race-specific rows (Asian,
    Black, Hispanic, White) had used counts explicitly noted as
    "Calculated N and X from sum of cases and controls" — i.e. cases
    plus matched controls summed together, giving exactly double the
    correct N. The scoping review's own Table 1 describes this study's
    population as "matched controls" only. Replaced with the
    controls-only counts and recomputed 95% CIs with the same
    Clopper-Pearson/Beta method used elsewhere in the dataset.
  - **ButCanc2020**: both race-specific rows (Black, White) are flagged
    ineligible by the co-author ("no indication of ethnicity" — SEER
    race/ethnicity was never confirmed). Dropped.
  - **VarCanc2020**: the "White" row was never confirmed as
    Non-Hispanic White (its SEER field was left as "Combined") and is
    flagged ineligible. Dropped; the (correctly flagged) "Black" row is
    unaffected.

  After this fix, the analytic sample produced by the existing
  filtering pipeline in `all_models_cluster.R` /
  `all_models_sensitivity_analysis_cluster.R` matches the Excel file's
  `viz_flag_race_eth_seer_comb == 1` rows exactly (previously 7 rows
  disagreed).

### Bug fixes

- **`model_code/all_models_sensitivity_analysis_cluster.R`**: fixed a
  bug where `n_chains` was referenced before it was ever assigned,
  causing the script to error out immediately (`object 'n_chains' not
  found`) before any model could be fit.
- **`model_code/LOO_all_model_comparison.R`**: fixed `fit_*.rds` paths
  that pointed at a `./cluster_output/` directory that neither fitting
  script writes to (they use `./results_cluster/` and
  `./results_cluster_sensitivity_analysis/`); fixed a `study_id` label
  inconsistency (`"NomJID2002"` vs `"NomJInf2002"` used everywhere
  else) that could silently misalign fitted study-level random effects
  during LOO recomputation; removed a leftover per-iteration debug
  `print()`.
- **`figures_code/figures_cluster_sep4.R`**: fixed a crashing bug where
  the max-prevalence risk-ratio/difference summary functions indexed a
  `"NH White"` matrix row that didn't exist (rows are named `"NHW"`),
  which errored out on every model; added the missing
  `figures_aug20_2026/binary_cohort_figs/` output directory creation;
  fixed the same `./cluster_output/` path issue; fixed hardcoded
  MCMC-chain/iteration counts in the trace-plot helper that silently
  mislabeled every trace plot's x-axis.
- **`figures_code/figures_cluster_sensitivity_analysis_sep15.R`**:
  removed a stray, copy-pasted block at the start of the "Model 1"
  section that referenced results objects before they existed and
  would have errored on a fresh run; fixed a filename bug where the
  "AI/AN" race label's `/` was passed through into a file path (via an
  incomplete `gsub`), which would break `ggsave()` for that group's
  per-race diagnostic plots; same `./cluster_output/` and trace-plot
  fixes as above. `cohort_years_to_plot` was used in the "Model 3 (RW
  cohort, expanded)" section but never defined, causing that section
  to error; per author decision, cohort-year subsetting was removed
  entirely from that section's max-prevalence and median-age-of-
  infection plots rather than restored, so they now show all cohort
  years.

None of these fixes change the JAGS model definitions or equations
themselves — all are either scripts that previously could not run to
completion, or presentation/robustness fixes.

### Author-driven changes made alongside this review

- **Force-of-infection intercept prior**: the Supplementary Appendix
  documents the γ-intercept prior as `t₍₃₎(−4.6, 1)` (Equations 9, 34,
  36), but the code used `dt(-3.2, 1, 3)` in JAGS across every model
  (binary and continuous cohort) — a real discrepancy between the
  published methods and the code that actually produced the results.
  Flagged during review; the author subsequently set the prior center
  to `-3.9` in `all_models_cluster.R`, which has been propagated to all
  three models in `all_models_sensitivity_analysis_cluster.R` so the
  same prior is used consistently across binary- and continuous-cohort
  formulations, per the Appendix's stated intent ("We assume similar
  priors for the alpha and gamma intercepts as for the corresponding
  fixed effects in the continuous and binary cohort base models"). The
  Supplementary Appendix text (which still says `-4.6`) should be
  updated to `-3.9` to match before the next submission/revision.

### Other findings (no code change)

- **`race_ethnicity_notes` / smaller row-level discrepancies**: the
  systematic Excel-vs-`.Rdata` comparison above covered every study;
  no other row-level mismatches were found beyond the three studies
  listed. `MuhJInf2014` is correctly absent from the modeling dataset
  (co-author-flagged ineligible, "Combined" only); the manually
  appended `NomJInf2002` row (patched in via `second_Nom_study` in all
  three model/figures scripts, since it's missing from the base
  `.Rdata`) was verified to exactly match the Excel record.
- **Force-of-infection at age 0**: `run_predictions()` in the
  continuous-cohort figures script sets `foi_mat[1, ] <- NA` (excludes
  age 0 from FOI plots), while the equivalent line in the binary-cohort
  figures script is commented out (age 0 is included). This is a minor
  inconsistency in presentation between the two figure sets, not a
  correctness bug.
