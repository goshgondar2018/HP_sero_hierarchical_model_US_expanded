# All code and data needed to run analyses conducted in "Widening racial and ethnic disparities in Helicobacter pylori acquisition in the USA, 1965–2014: a Bayesian meta-analytic catalytic model of age-specific force of infection"

See [NEWS.md](NEWS.md) for a changelog of fixes made to this repository since the original upload, including a data-correctness pass against a co-author-curated dataset and a number of crashing bugs in the figures/LOO scripts.

## 1. data

- df_hp_usa_w_final.Rdata: Meta-analytic serological dataset drawing from data provided in McMahon et al., 2025 and select additional studies. Dataset includes multiple fields related to study-specific features & results. HP denotes # testing positive for each subgroup_n by study.
- HPintheUS_WorkingFile_V12_20230727.xlsx: Co-author-curated working file behind the scoping review (McMahon et al., 2025), used to validate `df_hp_usa_w_final.Rdata`. The "Long HP data" sheet's `viz_flag_race_eth_seer_comb` column is the inclusion flag for the five-group SEER race/ethnicity categorization used in this repository's models (1 = include, 0 = exclude).

## 2. model_code

- all_models_cluster.R: Script to run all hierarchical models for binary cohort analysis (base model ("base"), base model with common kappa ("kappa_global"), base model with group-specific kappas ("kappa_partial"), base model with study random effects ("study_re"), base model with study random effects and common kappa ("kappa_global_study_re"), base model with study random effects and group-specific kappas ("kappa_partial_study_re")). Run `Rscript all_models_cluster.R [model_name ...]` from the repo root; with no arguments, all six models are run. Results (fit objects, convergence summaries, model comparison table) are written to `./results_cluster/`.
- all_models_sensitivity_analysis_cluster.R: Script to run all hierarchical models for continuous cohort analysis (base model ("base_continuous_bc"), random walk model with group-specific intercepts and common trajectories ("RW_continuous_bc"), random walk model with group-specific intercepts and trajectories ("RW_expanded_continuous_bc")). Same usage pattern as above; results are written to `./results_cluster_sensitivity_analysis/`.
- LOO_all_model_comparison.R: Script to generate out-of-sample model comparison metrics (PSIS-LOO) for all binary and continuous cohort models. Must be run after both scripts above, since it reads their saved fit objects from `./results_cluster/` and `./results_cluster_sensitivity_analysis/`.

## 3. figures_code

- figures_cluster_sep4.R: Script to run all figures and summary tables for the binary cohort analysis (reads fit objects from `./results_cluster/`, writes to `./figures_aug20_2026/binary_cohort_figs/` and `./results/`).
- figures_cluster_sensitivity_analysis_sep15.R: Script to run all figures and summary tables for the continuous cohort analysis, covering all three continuous-cohort models in one script (reads fit objects from `./results_cluster_sensitivity_analysis/`, writes to `./figures_aug20_2026/continuous_cohort_figs/` and `./results/`).

Both figures scripts must be run after the corresponding model_code script(s) above have produced fit objects.
