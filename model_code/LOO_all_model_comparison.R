

# *****************************************************************************
#### Load necessary packages ####
# *****************************************************************************

library(loo)
library(extraDistr)   
library(tidyverse)
# Required: supplies as.matrix()/varnames() methods for coda mcmc.list objects.
# Without it, as.matrix(mcmc_list) silently falls through to as.matrix.default
# and returns an unnamed list-matrix, so every subsequent
# mat_all[, "par[1,1]"] lookup fails with "no 'dimnames' attribute for array".
library(coda)

ilogit <- function(x) 1 / (1 + exp(-x))
race_labels <- c("Hispanic", "NH American Indian", "NH Asian", "NH Black", "NH White")
n_races <- length(race_labels)


dir.create("results", showWarnings = FALSE)

# *****************************************************************************
#### Indput data prep ####
# *****************************************************************************
load("data/df_hp_usa_w_final_unweighted_nhanes.Rdata")

second_Nom_study=data.frame(study_id="NomJInf2002",study_n=261.00,subgroup="Japanese",subgroup_race_ethnicity="Japanese",
subgroup_n=261.00,age=72.5,time_period=1982,birth_cohort=1909.5,seer_race_ethnicity="Non-Hispanic Asian or Pacific Islander",
race_ethnicity_notes=NA,hp_num=174+19,hp_denom=261.00,hp_pct=73.9,hp_notes=NA,version=NA,immigrant_n=0.16*261,year_pub=2002,
year=NA,race_f="Non-Hispanic Asian",total=261.00,HP=174+19,prev=74,prev_ci_lb=68,prev_ci_ub=79,birthplaceUS=NA)

data_model <- df_hp_usa_w_final %>%
  arrange(race_f, year) %>%
  filter(race_f != "Other") 

data_model_complete=rbind.data.frame(data_model,second_Nom_study)%>%
  mutate(
    study_id_fac = factor(study_id),
    race_fac     = factor(race_f),
    time         = year - 1988,
    bc_scaled    = birth_cohort - min(birth_cohort),
    birth_cohort_cat = ifelse(birth_cohort < 1960, 0, 1)
  )

data_model_filtered <- data_model_complete %>%
  filter(study_id != 'CorCanc1990') %>%
  filter(study_id != 'GraGast1991') %>%
  filter(study_id != 'SawJRhe2004') %>%
  filter(study_id != 'MuhJInf2015') %>%
  filter(study_id != 'EppCanc2011_3') %>%
  filter(!(study_id == 'GooJImm2003' & subgroup_n == 384.00)) %>%
  filter(!(study_id == 'FelBirt2007' & subgroup_n == 147.00)) %>%
  filter(!(study_id == 'ParClin2000' & subgroup_n == 2080.00)) %>%
  filter(!(study_id == 'ParClin2000' & subgroup_n == 1040.00))

N          <- nrow(data_model_filtered)
inds       <- as.numeric(data_model_filtered$race_fac)
n_races    <- length(unique(data_model_filtered$race_fac))
hp         <- round(data_model_filtered$HP)
age        <- data_model_filtered$age
bc_scaled  <- data_model_filtered$bc_scaled
birth_cohort_cat <- data_model_filtered$birth_cohort_cat
# Binomial denominator: the number ASSAYED, not the number enrolled.
#
# `subgroup_n` is the subgroup size; `total` (== `hp_denom`) is how many were
# actually tested. They differ for three rows, and the dataset's own `prev`
# column is computed from `hp_denom`, so using subgroup_n made the modelled
# prevalence disagree with the reported prevalence for those rows:
#
#   MeiSocS2016 Hispanic  1149/1562 = 73.6%  vs  1149/1226 = 93.7%  (20.2 pp)
#   RubGenE2013 AI/AN      594/801  = 74.2%  vs   594/766  = 77.5%  ( 3.4 pp)
#   MieHeli2018 AI/AN      422/571  = 73.9%  vs   422/561  = 75.2%  ( 1.3 pp)
#
# MeiSocS2016 is a 1,226-person Hispanic observation in the pre-1960 cohort, so
# the old denominator biased Hispanic maximum prevalence downward.
n          <- round(data_model_filtered$total)
study_id_raw <- data_model_filtered$study_id
study_id     <- as.numeric(as.factor(study_id_raw))
n_studies    <- length(unique(study_id))
cohort_idx <- as.numeric(as.factor(round(data_model_filtered$birth_cohort)))
n_cohorts <- length(unique(cohort_idx))

# ******************************************************************************************
#### Compute p[i] matrix [n_iter x N] for a given model structure/kappa/study_re combo ####
# ******************************************************************************************

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
    p_mat[it, ] <- 0.001 + 0.998 * p_raw
  }

  list(p_mat = p_mat, theta_draws = theta_draws)
}

# ******************************************************************************
#### Analytic log-likelihood ####
# ******************************************************************************

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

# ******************************************************************************
# Run LOO for one model and return log-lik matrix
# ******************************************************************************
run_loo_for_model <- function(model_name, fit_path, p_mat_args, hp, n) {

  fit_obj <- readRDS(fit_path)
  mcmc_list <- fit_obj$mcmc_all         
  mat_all <- as.matrix(mcmc_list)

  pm <- do.call(compute_p_matrix, c(list(mat_all = mat_all), p_mat_args))
  loglik <- compute_loglik_betabinom(pm$p_mat, pm$theta_draws, hp, n)

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

# *****************************************************************************
# Registry of all models
# *****************************************************************************
loo_configs <- list(
  base = list(
    file = "./results_cluster/fit_base.rds",
    args = list(structure = "linear_binary", kappa_mode = "none",
                has_study_re = FALSE, bc = birth_cohort_cat,
                inds = inds, age = age, n_races = n_races)
  ),
  study_re = list(
    file = "./results_cluster/fit_study_re.rds",
    args = list(structure = "linear_binary", kappa_mode = "none",
                has_study_re = TRUE, study_id = study_id, bc = birth_cohort_cat,
                inds = inds, age = age, n_races = n_races)
  ),
  kappa_global = list(
    file = "./results_cluster/fit_kappa_global.rds",
    args = list(structure = "linear_binary", kappa_mode = "global",
                has_study_re = FALSE, bc = birth_cohort_cat,
                inds = inds, age = age, n_races = n_races)
  ),
  kappa_partial = list(
    file = "./results_cluster/fit_kappa_partial.rds",
    args = list(structure = "linear_binary", kappa_mode = "partial",
                has_study_re = FALSE, bc = birth_cohort_cat,
                inds = inds, age = age, n_races = n_races)
  ),
  kappa_global_study_re = list(
    file = "./results_cluster/fit_kappa_global_study_re.rds",
    args = list(structure = "linear_binary", kappa_mode = "global",
                has_study_re = TRUE, study_id = study_id, bc = birth_cohort_cat,
                inds = inds, age = age, n_races = n_races)
  ),
  kappa_partial_study_re = list(
    file = "./results_cluster/fit_kappa_partial_study_re.rds",
    args = list(structure = "linear_binary", kappa_mode = "partial",
                has_study_re = TRUE, study_id = study_id, bc = birth_cohort_cat,
                inds = inds, age = age, n_races = n_races)
  ),
  base_continuous_bc = list(
    file = "./results_cluster_sensitivity_analysis/fit_base_continuous_bc.rds",
    args = list(structure = "linear_continuous", kappa_mode = "none",
               has_study_re = FALSE, bc_scaled = bc_scaled,
               inds = inds, age = age, n_races = n_races)
 ),
  RW_continuous_bc = list(
    file = "./results_cluster_sensitivity_analysis/fit_RW_continuous_bc.rds",
  args = list(structure = "rw_shared", kappa_mode = "none",
                has_study_re = TRUE, study_id = study_id,
              cohort_idx = cohort_idx, n_cohorts = n_cohorts,
               inds = inds, age = age, n_races = n_races)
  ),
  RW_expanded_continuous_bc = list(
    file = "./results_cluster_sensitivity_analysis/fit_RW_expanded_continuous_bc.rds",
    args = list(structure = "rw_race_dev", kappa_mode = "none",
               has_study_re = TRUE, study_id = study_id,
               cohort_idx = cohort_idx, n_cohorts = n_cohorts,
              inds = inds, age = age, n_races = n_races)
  )
)


# *****************************************************************************
# Run all models
# *****************************************************************************
models_to_run <- names(loo_configs)   

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

# *****************************************************************************
# Generate model comparison table 
# *****************************************************************************
loo_list <- lapply(loo_results, function(x) x$loo_result)
comparison <- loo::loo_compare(loo_list)

comparison_df <- as.data.frame(comparison) %>%
  rownames_to_column("model") %>%
  arrange(desc(elpd_loo))

write_csv(comparison_df, "results/loo_comparison_ALL.csv")

fit_obj <- readRDS("./results_cluster/fit_study_re.rds")
mat_all <- as.matrix(fit_obj$mcmc_all)
study_re_cols <- grep("^study_re\\[", colnames(mat_all))
cat("n_studies (in saved fit):", length(study_re_cols), "\n")