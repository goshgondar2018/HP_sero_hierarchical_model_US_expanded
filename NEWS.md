# NEWS

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
  would have errored on a fresh run; defined `cohort_years_to_plot`,
  which was used in the "Model 3 (RW cohort, expanded)" section but
  never defined, causing that section to error; fixed a filename bug
  where the "AI/AN" race label's `/` was passed through into a file
  path (via an incomplete `gsub`), which would break `ggsave()` for
  that group's per-race diagnostic plots; same `./cluster_output/` and
  trace-plot fixes as above.

None of these fixes change the JAGS model definitions, priors, or
equations themselves — all are either scripts that previously could not
run to completion, or presentation/robustness fixes.

### Known issues flagged for author review (not changed)

These were surfaced during the review but intentionally left
unchanged, since resolving them one way or the other affects reported
results and should be a deliberate author decision:

- **Force-of-infection intercept prior mismatch**: the Supplementary
  Appendix documents the γ-intercept prior as `t₍₃₎(−4.6, 1)` (Equations
  9, 34, 36), but every model script (binary and continuous cohort) uses
  `dt(-3.2, 1, 3)` in JAGS. This is a real discrepancy between the
  published methods and the code that actually produced the results —
  worth confirming which is correct and updating the other to match.
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
