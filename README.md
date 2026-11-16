# All code and data needed to run analyses conducted in "Widening racial and ethnic disparities in Helicobacter pylori acquisition in the USA, 1965–2014: a Bayesian meta-analytic catalytic model of age-specific force of infection"

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


## 1. data

- df_hp_usa_w_final_unweighted_nhanes.Rdata: Meta-analytic serological dataset drawing from data provided in McMahon et al., 2025 and select additional studies. Dataset includes multiple fields related to study-specific features & results. HP denotes # testing positive by study subgroup. `total` (== `hp_denom`) is the number actually assayed and is the binomial denominator used by the models. Use unweighted counts for NHANES studies
- HPintheUS_WorkingFile_V12_20230727.xlsx: Co-author-curated working file behind the scoping review (McMahon et al., 2025), used to validate `df_hp_usa_w_final.Rdata`. The "Long HP data" sheet's `viz_flag_race_eth_seer_comb` column is the inclusion flag for the five-group SEER race/ethnicity categorization used in this repository's models (1 = include, 0 = exclude)

## 2. model_code

- all_models_cluster.R: Script to run all hierarchical models for binary cohort analysis (base model ("base"), base model with common kappa ("kappa_global"), base model with group-specific kappas ("kappa_partial"), base model with study random effects ("study_re"), base model with study random effects and common kappa ("kappa_global_study_re"), base model with study random effects and group-specific kappas ("kappa_partial_study_re")). Running the script with no arguments runs all six models
- all_models_sensitivity_analysis_cluster.R: Script to run all hierarchical models for continuous cohort analysis (base model ("base_continuous_bc"), random walk model with group-specific intercepts and common trajectories ("RW_continuous_bc"), random walk model with group-specific intercepts and trajectories ("RW_expanded_continuous_bc")). Running the script with no arguments runs all three models
- LOO_all_model_comparison.R: Script to generate out-of-sample model comparison metrics (PSIS-LOO) for all binary and continuous cohort models. Must be run after both scripts above, since it reads their saved fit objects

## 3. figures_code

- figures_cluster_sep4.R: Script to run all figures and summary tables for the binary cohort analysis
- figures_cluster_sensitivity_analysis_sep15.R: Script to run all figures and summary tables for the continuous cohort analysis

