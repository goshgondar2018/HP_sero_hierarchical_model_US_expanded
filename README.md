# All code and data needed to run analyses conducted in "Widening racial and ethnic disparities in Helicobacter pylori acquisition in the USA, 1965–2014: a Bayesian meta-analytic catalytic model of age-specific force of infection"

## 1. data

df_hp_usa_w_final.Rdata: Meta-analytic serological dataset drawing from data provided in McMahon et al., 2025 and select additional studies. Dataset includes multiple fields related to study-specific features & results. HP denotes # testing positive for each subgroup_n by study. 

## 2. model_code

- all_models_cluster.R: Script to run all hierarchical models for binary cohort analysis (base model ("base"), base model with common kappa ("kappa_global"), base model with group-specific kappas ("kappa_partial"), base model with study random effects ("study_re"), base model with study random effects and group-specific kappas ("kappa_partial_study_re"))
- all_models_sensitivity_analysis_cluster.R: Script to run all hierarchical models for continuous cohort analysis (base model ("base_continuous_bc", random walk model with group-specific intercepts and common trajectories ("RW_continuous_bc"), random walk model with group-specific intercepts and trajectories ("RW_expanded_continuous_bc"))
- LOO_all_model_comparison.R: Script to generate out-of-sample model comparison metrics for all binary and continuous cohort models

## 3. figures_code

- figures_cluster_sep4.R: Script to run all figures summarizing binary cohort analysis results
- figures_cluster_sensitivity_analysis_base_model.R: Script to run all figures summarizing continuous cohort base model analysis results
- figures_cluster_sensitivity_analysis_RW1.R: Script to run all figures summarizing continuous cohort RW_continuous_bc model analysis results
- figures_cluster_sensitivity_analysis_RWexpanded.R: Script to run all figures summarizing continuous cohort RW_expanded_continuous_bc model analysis results
