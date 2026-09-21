# All code and data needed to run analyses conducted in "Widening racial and ethnic disparities in Helicobacter pylori acquisition in the USA, 1965–2014: a Bayesian meta-analytic catalytic model of age-specific force of infection"

> **New and improved version.** The data and code have been reviewed twice
> against the co-author-curated working file and the Supplementary Appendix.
> Several corrections change reported results — most importantly, the force of
> infection was previously understated by up to 27% at the ages where it
> peaks, which shifts both peak FOI and the age at peak, differently for each
> race/ethnicity group. **Read [NEWS.md](NEWS.md) before using these results or
> comparing them with the submitted manuscript.** All models need refitting and
> all results regenerating.

## Requirements

- **R** ≥ 4.1 (developed against 4.5.1)
- **JAGS** ≥ 4.3, installed separately from R
  (`brew install jags`, `apt install jags`, or
  [mcmc-jags.sourceforge.io](https://mcmc-jags.sourceforge.io))
- R packages:

```r
install.packages(c("tidyverse", "data.table", "runjags", "coda", "rjags",
                   "matrixStats", "gt", "loo", "extraDistr", "reshape2"))
```

Optional: `webshot2` + `chromote` render the PSRF diagnostic tables as PNG.
Without them those tables are still written as CSV.

**Run every script from the repository root** — all paths are relative to it:

```sh
cd /path/to/HP_sero_hierarchical_model_US_expanded
Rscript model_code/all_models_cluster.R
```

## Order of operations

```sh
# 1. Fit the models (hours -- see runtime note below)
Rscript model_code/all_models_cluster.R                     # all six binary-cohort models
Rscript model_code/all_models_sensitivity_analysis_cluster.R # all three continuous-cohort models

# 2. Out-of-sample model comparison (needs step 1)
Rscript model_code/LOO_all_model_comparison.R

# 3. Figures and summary tables (needs step 1)
Rscript figures_code/figures_cluster_sep4.R
Rscript figures_code/figures_cluster_sensitivity_analysis_sep15.R
```

Each script also accepts model names to fit or plot a subset, e.g.
`Rscript model_code/all_models_cluster.R base study_re`.

**Runtime.** 3 chains × (10,000 adaptation + 500,000 burn-in + 20,000 retained
samples, thinned by 20). Each model takes hours; chains run in parallel, one
per core.

## 1. data

- df_hp_usa_w_final.Rdata: Meta-analytic serological dataset drawing from data provided in McMahon et al., 2025 and select additional studies. Dataset includes multiple fields related to study-specific features & results. HP denotes # testing positive by study subgroup. `total` (== `hp_denom`) is the number actually assayed and is the binomial denominator used by the models; `subgroup_n` is the number enrolled and differs from it for three rows (see [NEWS.md](NEWS.md))
- HPintheUS_WorkingFile_V12_20230727.xlsx: Co-author-curated working file behind the scoping review (McMahon et al., 2025), used to validate `df_hp_usa_w_final.Rdata`. The "Long HP data" sheet's `viz_flag_race_eth_seer_comb` column is the inclusion flag for the five-group SEER race/ethnicity categorization used in this repository's models (1 = include, 0 = exclude)

## 2. model_code

- all_models_cluster.R: Script to run all hierarchical models for binary cohort analysis (base model ("base"), base model with common kappa ("kappa_global"), base model with group-specific kappas ("kappa_partial"), base model with study random effects ("study_re"), base model with study random effects and common kappa ("kappa_global_study_re"), base model with study random effects and group-specific kappas ("kappa_partial_study_re")). Running the script with no arguments runs all six models
- all_models_sensitivity_analysis_cluster.R: Script to run all hierarchical models for continuous cohort analysis (base model ("base_continuous_bc"), random walk model with group-specific intercepts and common trajectories ("RW_continuous_bc"), random walk model with group-specific intercepts and trajectories ("RW_expanded_continuous_bc")). Running the script with no arguments runs all three models
- LOO_all_model_comparison.R: Script to generate out-of-sample model comparison metrics (PSIS-LOO) for all binary and continuous cohort models. Must be run after both scripts above, since it reads their saved fit objects

## 3. figures_code

- figures_cluster_sep4.R: Script to run all figures and summary tables for the binary cohort analysis
- figures_cluster_sensitivity_analysis_sep15.R: Script to run all figures and summary tables for the continuous cohort analysis

## Outputs

| Location | Contents |
|---|---|
| `results_cluster/` | Fitted binary-cohort models (`fit_<model>.rds`), per-model summaries, `model_comparison_ALL.csv` |
| `results_cluster_sensitivity_analysis/` | Fitted continuous-cohort models and their comparison table |
| `results/` | Prevalence and FOI curves, max prevalence, prevalence ratios and differences, index of disparity, median age of infection, PSIS-LOO comparison |
| `figures_aug20_2026/binary_cohort_figs/` | Binary-cohort figures and PSRF tables |
| `figures_aug20_2026/continuous_cohort_figs/` | Continuous-cohort figures and PSRF tables |

These directories are git-ignored; regenerate them by running the pipeline.

## Checking convergence before reporting

Each fitted object records `mean_gelman_all`, `max_gelman_all`,
`n_PSRF_above_1.1` and a `converged` flag, and the per-model summary CSVs carry
the same columns. Judge convergence on the **maximum** PSRF, not the mean — a
single badly mixing parameter is invisible in an average over 15–40
parameters. The fitting scripts print `*** NOT CONVERGED ***` and raise an R
warning when the maximum exceeds 1.1.

DIC and ELPD-LOO are not on the same scale: DIC comes from the JAGS `deviance`
node, which conditions on the latent `pi[i]`, while the PSIS-LOO
log-likelihood integrates `pi[i]` out analytically. ELPD-LOO is the
appropriate basis for model comparison. Both are meaningless if the chains
have not mixed.

## Prior parameterisation

JAGS `dnorm(mu, tau)` and `dt(mu, tau, k)` take a **precision**, not a standard
deviation or a variance. Three priors were originally written as though the
second argument were an SD or a variance. The γ-intercept has been resolved to
`-3.9` by author decision. The other two are exposed as named constants so
either reading can be run without editing model code:

| Constant | File | Default | Alternative |
|---|---|---|---|
| `KAPPA_TAU` | `all_models_cluster.R` | `0.40` (SD 1.58, as coded) | `6.25` (SD 0.40, Appendix reading) |
| `RW_INIT_TAU` | `all_models_sensitivity_analysis_cluster.R` | `25.00` (SD 0.20, as coded) | `0.04` (SD 5, Appendix variance 25) |

Both are documented inline with the implied prior intervals. See
[NEWS.md](NEWS.md) for why each matters.
