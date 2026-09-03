# =============================================================================
# PSIS-LOO (approximate leave-one-out cross-validation) for all fitted models
#
# Unlike the PPC/RMSE/CRPS checks (in-sample: model checked against data it
# was fit on), PSIS-LOO approximates genuine leave-one-out predictive
# performance -- "how well would this model predict each observation if it
# hadn't seen it" -- computed from the SAME posterior draws, no refitting.
#
# Requires pointwise log-likelihood per observation, per draw. Since the
# beta-binomial likelihood has a closed-form marginal density (integrating
# out the latent pi[i] analytically), we compute this directly rather than
# needing pi[i] monitored in JAGS:
#
#   y[i] | p[i], theta  ~ BetaBinomial(n[i], a[i]=p[i]*theta, b[i]=(1-p[i])*theta)
#
#   log BetaBinom(y; n, a, b) = log C(n,y) + lbeta(y+a, n-y+b) - lbeta(a,b)
# =============================================================================

if (!requireNamespace("loo", quietly = TRUE)) install.packages("loo")
if (!requireNamespace("extraDistr", quietly = TRUE)) install.packages("extraDistr")
library(loo)
library(extraDistr)   # provides dbbinom() -- beta-binomial density
library(tidyverse)

ilogit <- function(x) 1 / (1 + exp(-x))
race_labels <- c("Hispanic", "NH American Indian", "NH Asian", "NH Black", "NH White")
n_races <- length(race_labels)

if (!exists("cohort_idx")) {
  cohort_idx <- as.numeric(as.factor(round(data_model_filtered$birth_cohort)))
}
n_cohorts <- length(unique(cohort_idx))

dir.create("results", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Compute p[i] matrix [n_iter x N] for a given model structure/kappa/study_re
# combo. Mirrors the simulate_yrep logic exactly, but returns p[i] draws
# directly (no simulation) -- reused for the analytic log-likelihood.
# -----------------------------------------------------------------------------
compute_p_matrix <- function(mat_all, age, inds,
                              bc = NULL, bc_scaled = NULL, cohort_idx = NULL,
                              structure = c("linear_binary", "linear_continuous",
                                            "rw_shared", "rw_race_dev"),
                              kappa_mode = c("none", "global", "partial"),
                              has_study_re = FALSE, study_id = NULL,
                              n_races, n_cohorts = NULL) {

  structure  <- match.arg(structure)
  kappa_mode <- match.arg(kappa_mode)
  N <- length(age)
  n_iter <- nrow(mat_all)

  study_re_draws <- NULL
  if (has_study_re) {
    n_studies <- length(unique(study_id))
    study_re_draws <- sapply(1:n_studies, function(s) mat_all[, paste0("study_re[", s, "]")])
  }

  m_kappa <- NULL
  if (kappa_mode == "global") {
    m_kappa <- matrix(mat_all[, "kappa"], nrow = n_iter, ncol = n_races)
  } else if (kappa_mode == "partial") {
    m_kappa <- mat_all[, paste0("kappa[", 1:n_races, "]"), drop = FALSE]
  }

  if (structure %in% c("linear_binary", "linear_continuous")) {
    get_par_matrix <- function(param_col) {
      mat_all[, paste0("par[", 1:n_races, ",", param_col, "]"), drop = FALSE]
    }
    m_alpha0 <- get_par_matrix(1); m_alpha1 <- get_par_matrix(2)
    m_gamma0 <- get_par_matrix(3); m_gamma1 <- get_par_matrix(4)
  } else {
    m_alpha0 <- mat_all[, paste0("par_alpha0[", 1:n_races, "]"), drop = FALSE]
    m_gamma0 <- mat_all[, paste0("par_gamma0[", 1:n_races, "]"), drop = FALSE]
    rw_alpha_pop <- mat_all[, paste0("rw_alpha_pop[", 1:n_cohorts, "]"), drop = FALSE]
    rw_gamma_pop <- mat_all[, paste0("rw_gamma_pop[", 1:n_cohorts, "]"), drop = FALSE]
    if (structure == "rw_race_dev") {
      rw_alpha_dev <- array(NA, dim = c(n_iter, n_races, n_cohorts))
      rw_gamma_dev <- array(NA, dim = c(n_iter, n_races, n_cohorts))
      for (j in 1:n_races) for (t in 1:n_cohorts) {
        rw_alpha_dev[, j, t] <- mat_all[, paste0("rw_alpha_dev[", j, ",", t, "]")]
        rw_gamma_dev[, j, t] <- mat_all[, paste0("rw_gamma_dev[", j, ",", t, "]")]
      }
    }
  }

  theta_draws <- mat_all[, "theta"]
  p_mat <- matrix(NA, nrow = n_iter, ncol = N)

  for (it in 1:n_iter) {
    if (structure == "linear_binary") {
      a_logit <- m_alpha0[it, inds] + m_alpha1[it, inds] * bc
      g_log   <- m_gamma0[it, inds] + m_gamma1[it, inds] * bc
    } else if (structure == "linear_continuous") {
      a_logit <- m_alpha0[it, inds] + m_alpha1[it, inds] * bc_scaled
      g_log   <- m_gamma0[it, inds] + m_gamma1[it, inds] * bc_scaled
    } else if (structure == "rw_shared") {
      a_logit <- m_alpha0[it, inds] + rw_alpha_pop[it, cohort_idx]
      g_log   <- m_gamma0[it, inds] + rw_gamma_pop[it, cohort_idx]
    } else {
      dev_a <- rw_alpha_dev[it, , ][cbind(inds, cohort_idx)]
      dev_g <- rw_gamma_dev[it, , ][cbind(inds, cohort_idx)]
      a_logit <- m_alpha0[it, inds] + rw_alpha_pop[it, cohort_idx] + dev_a
      g_log   <- m_gamma0[it, inds] + rw_gamma_pop[it, cohort_idx] + dev_g
    }
    if (has_study_re) a_logit <- a_logit + study_re_draws[it, study_id]

    alpha_i <- ilogit(a_logit)
    gamma_i <- exp(g_log)

    if (kappa_mode == "none") {
      p_raw <- alpha_i * (1 - exp(-gamma_i * age))
    } else {
      kappa_i <- m_kappa[it, inds]
      p_raw <- alpha_i * (1 - exp(-(gamma_i * age)^kappa_i))
    }
    print(it)
    p_mat[it, ] <- 0.001 + 0.998 * p_raw
  }

  list(p_mat = p_mat, theta_draws = theta_draws)
}

# -----------------------------------------------------------------------------
# Analytic pointwise log-likelihood for the beta-binomial marginal
# -----------------------------------------------------------------------------
compute_loglik_betabinom <- function(p_mat, theta_draws, hp, n) {
  n_iter <- nrow(p_mat)
  N <- ncol(p_mat)
  loglik <- matrix(NA, nrow = n_iter, ncol = N)
  for (it in 1:n_iter) {
    a <- p_mat[it, ] * theta_draws[it]
    b <- (1 - p_mat[it, ]) * theta_draws[it]
    loglik[it, ] <- extraDistr::dbbinom(hp, size = n, alpha = a, beta = b, log = TRUE)
  }
  loglik
}

# -----------------------------------------------------------------------------
# Run LOO for one model: builds log-lik matrix, computes relative_eff from
# chain structure (for more reliable PSIS diagnostics), runs loo()
# -----------------------------------------------------------------------------
run_loo_for_model <- function(model_name, fit_path, p_mat_args, hp, n) {
  cat("Computing LOO for:", model_name, "\n")

  fit_obj <- readRDS(fit_path)
  mcmc_list <- fit_obj$mcmc_all         # coda mcmc.list -- preserves chain structure
  mat_all <- as.matrix(mcmc_list)

  pm <- do.call(compute_p_matrix, c(list(mat_all = mat_all), p_mat_args))
  loglik <- compute_loglik_betabinom(pm$p_mat, pm$theta_draws, hp, n)

  # chain_id for relative_eff -- matches as.matrix.mcmc.list's row order
  # (chains stacked sequentially)
  n_chains <- length(mcmc_list)
  n_iter_chain <- nrow(mcmc_list[[1]])
  chain_id <- rep(1:n_chains, each = n_iter_chain)

  r_eff <- tryCatch(
    loo::relative_eff(exp(loglik), chain_id = chain_id, cores = 1),
    error = function(e) { cat("  relative_eff failed, proceeding without it:", conditionMessage(e), "\n"); NULL }
  )

  loo_result <- loo::loo(loglik, r_eff = r_eff)
  list(model_name = model_name, loo_result = loo_result, loglik = loglik)
}

# -----------------------------------------------------------------------------
# Registry: all models, both families, with their p_mat_args
# -----------------------------------------------------------------------------
loo_configs <- list(
  base = list(
    file = "./cluster_output/fit_base.rds",
    args = list(structure = "linear_binary", kappa_mode = "none",
                has_study_re = FALSE, bc = birth_cohort_cat,
                inds = inds, age = age, n_races = n_races)
  ),
  study_re = list(
    file = "./cluster_output/fit_study_re.rds",
    args = list(structure = "linear_binary", kappa_mode = "none",
                has_study_re = TRUE, study_id = study_id, bc = birth_cohort_cat,
                inds = inds, age = age, n_races = n_races)
  ),
  kappa_global = list(
    file = "./cluster_output/fit_kappa_global.rds",
    args = list(structure = "linear_binary", kappa_mode = "global",
                has_study_re = FALSE, bc = birth_cohort_cat,
                inds = inds, age = age, n_races = n_races)
  ),
  kappa_partial = list(
    file = "./cluster_output/fit_kappa_partial.rds",
    args = list(structure = "linear_binary", kappa_mode = "partial",
                has_study_re = FALSE, bc = birth_cohort_cat,
                inds = inds, age = age, n_races = n_races)
  ),
  kappa_partial_study_re = list(
    file = "./cluster_output/fit_kappa_partial_study_re.rds",
    args = list(structure = "linear_binary", kappa_mode = "partial",
                has_study_re = TRUE, study_id = study_id, bc = birth_cohort_cat,
                inds = inds, age = age, n_races = n_races)
  )#,
  #base_continuous_bc = list(
   # file = "./cluster_output/fit_base_continuous_bc.rds",
   # args = list(structure = "linear_continuous", kappa_mode = "none",
             #   has_study_re = FALSE, bc_scaled = bc_scaled,
            #    inds = inds, age = age, n_races = n_races)
  #),
  #RW_continuous_bc = list(
    #file = "./cluster_output/fit_RW_continuous_bc.rds",
    #args = list(structure = "rw_shared", kappa_mode = "none",
               # has_study_re = TRUE, study_id = study_id,
               # cohort_idx = cohort_idx, n_cohorts = n_cohorts,
               # inds = inds, age = age, n_races = n_races)
  #),
  #RW_expanded_continuous_bc = list(
   # file = "./cluster_output/fit_RW_expanded_continuous_bc.rds",
   # args = list(structure = "rw_race_dev", kappa_mode = "none",
            #    has_study_re = TRUE, study_id = study_id,
            #    cohort_idx = cohort_idx, n_cohorts = n_cohorts,
            #    inds = inds, age = age, n_races = n_races)
 # )
)

# -----------------------------------------------------------------------------
# Run all models -- adjust this list to just the ones you have .rds files for
# -----------------------------------------------------------------------------
models_to_run <- names(loo_configs)   # or specify a subset, e.g. c("base","study_re",...)

loo_results <- list()
for (model_name in models_to_run) {
  cfg <- loo_configs[[model_name]]
  loo_results[[model_name]] <- run_loo_for_model(model_name, cfg$file, cfg$args, hp, n)
  print(loo_results[[model_name]]$loo_result)

  # Pareto k diagnostic check -- flags unreliable PSIS estimates
  n_bad_k <- sum(loo_results[[model_name]]$loo_result$diagnostics$pareto_k > 0.7)
  if (n_bad_k > 0) {
    cat("  WARNING:", n_bad_k, "observations have Pareto k > 0.7 -- LOO estimate",
        "may be unreliable for those points.\n")
  }
}

# -----------------------------------------------------------------------------
# Model comparison table (elpd_diff, se_diff -- standard loo_compare output)
# -----------------------------------------------------------------------------
loo_list <- lapply(loo_results, function(x) x$loo_result)
comparison <- loo::loo_compare(loo_list)

cat("\n===== PSIS-LOO model comparison (best model first) =====\n")
print(comparison)

comparison_df <- as.data.frame(comparison) %>%
  rownames_to_column("model") %>%
  arrange(desc(elpd_loo))

write_csv(comparison_df, "results/loo_comparison_ALL.csv")

cat("\nInterpretation:\n",
    "- elpd_diff / se_diff: models are listed best-first (highest elpd_loo).\n",
    "  If |elpd_diff| < ~2*se_diff for two models, the difference in\n",
    "  out-of-sample predictive performance is NOT clearly distinguishable.\n",
    "- looic: LOO information criterion (-2*elpd_loo), lower is better --\n",
    "  directly comparable to DIC/WAIC in spirit, but a genuine approximation\n",
    "  of leave-one-out predictive accuracy rather than an in-sample penalty.\n",
    "- Check the Pareto k warnings above; many k>0.7 values means PSIS-LOO\n",
    "  itself may be unreliable for that model (rare here given moderate N).\n")


# Pareto-k diagnostics, specifically for base/kappa_global/kappa_partial --
# confirm the ~70-point gap isn't driven by a handful of unreliable PSIS points
sum(loo_results[["base"]]$loo_result$diagnostics$pareto_k > 0.7)
sum(loo_results[["kappa_partial"]]$loo_result$diagnostics$pareto_k > 0.7)
sum(loo_results[["kappa_global"]]$loo_result$diagnostics$pareto_k > 0.7)