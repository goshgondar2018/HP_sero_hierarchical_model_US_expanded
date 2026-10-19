
# *****************************************************************************
#### Load necessary packages ####
# *****************************************************************************

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr) 
  library(stringr)
  library(runjags)
  library(coda)
})

set.seed(123)


# *****************************************************************************
#### Parallel execution settings ####
# *****************************************************************************
n_chains <- 3
n_cores_avail <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores_avail)) n_cores_avail <- parallel::detectCores()
n_cores_use <- min(n_chains, max(1, n_cores_avail))
cat(sprintf("Using %d cores for %d chains.\n", n_cores_use, n_chains))

out_dir <- "./results_cluster_sensitivity_analysis"
fig_dir <- "./figures_cluster_sensitivity_analysis"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# *****************************************************************************
#### Prior parameterisation ####
RW_INIT_TAU <- 25.00
# *****************************************************************************

# *****************************************************************************
#### Data prep (for all models) ####
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

## analytic sample size
total_pop = sum(data_model_filtered$total)

N          <- nrow(data_model_filtered)
inds       <- as.numeric(data_model_filtered$race_fac)
n_races    <- length(unique(data_model_filtered$race_fac))
hp         <- round(data_model_filtered$HP)
age        <- data_model_filtered$age
bc_scaled  <- data_model_filtered$bc_scaled

n          <- round(data_model_filtered$total)
OmegaB     <- diag(4)

study_id_raw <- data_model_filtered$study_id
study_id     <- as.numeric(as.factor(study_id_raw))
n_studies    <- length(unique(study_id))

study_lookup <- data.frame(
  study_id_num = study_id,
  study_abbrev = as.character(study_id_raw)
) %>% distinct() %>% arrange(study_id_num)
saveRDS(study_lookup, file.path(out_dir, "study_lookup.rds"))

cohort_idx <- as.numeric(as.factor(round(data_model_filtered$birth_cohort)))
n_cohorts  <- length(unique(cohort_idx))

base_data_continuous <- list(
  y = hp, x = age, bc = bc_scaled,   
  n = n, N = N, inds = inds, OmegaB = OmegaB, n_races = n_races,
  n_studies  = n_studies, study_id  = study_id
)
rw_data <- c(base_data_continuous, list(
  cohort_idx = cohort_idx, n_cohorts = n_cohorts,
  rw_init_tau = RW_INIT_TAU
))


# *****************************************************************************
#### Model definitions ####
# *****************************************************************************
model_base_continuous_bc <- paste0("model{
  for (i in 1:N){
    y[i]   ~ dbinom(pi[i], n[i])
    pi[i]  ~ dbeta(a[i], b[i])
    a[i]   <- p[i] * theta
    b[i]   <- (1 - p[i]) * theta
    
    ## bound p away from 0 so log() values not NaN
    p[i]   <- 0.001 + 0.998 * p_raw[i]
    p_raw[i] <- alpha[i] * (1 - exp(-gamma[i] * x[i]))
    
    logit(alpha[i]) <- par[inds[i],1] + par[inds[i],2]*bc[i] + study_re[study_id[i]]
    log(gamma[i])   <- par[inds[i],3] + par[inds[i],4]*bc[i]
  }
  
  ## overdispersion parameter
  log_theta ~ dnorm(3, 1)   
  theta     <- exp(log_theta)
  
  ## race-specific random effects
  for (j in 1:n_races){
    par[j,1:4] ~ dmnorm(mub[1:4], taub[1:4,1:4]) 
  }

  ## study level random effects
  for (s in 1:n_studies){
    study_re[s] ~ dnorm(0, tau_study)
  }
  tau_study ~ dgamma(2, 0.5)
  sd_study  <- 1/sqrt(tau_study)

  ## fixed effects
  mub[1] ~ dt(0, 1, 1)
  mub[2] ~ dt(0, 1, 1)
  mub[3] ~ dt(-3.9, 1, 3)
  mub[4] ~ dt(0, 1, 1)
  
  taub[1:4,1:4] ~ dwish(OmegaB[1:4,1:4], 5)
  sigma2_b[1:4,1:4] <- inverse(taub[,])

  for(k in 1:4){
    var_b[k] <- sigma2_b[k,k]
  }
  
  cov_b[1] <- sigma2_b[2,1]
  cov_b[2] <- sigma2_b[3,1]
  cov_b[3] <- sigma2_b[3,2]
  cov_b[4] <- sigma2_b[4,1]
  cov_b[5] <- sigma2_b[4,2]
  cov_b[6] <- sigma2_b[4,3]
  
  ForMCMCConv[1:4]  <- mub
  ForMCMCConv[5:8]  <- var_b
  ForMCMCConv[9:14] <- cov_b
  ForMCMCConv[15]   <- theta
  ForMCMCConv[16]   <- sd_study

}")

model_RW_continuous_bc <- paste0("model{

  for (i in 1:N){
    y[i]   ~ dbinom(pi[i], n[i])
    pi[i]  ~ dbeta(a[i], b[i])
    a[i]   <- p[i] * theta
    b[i]   <- (1 - p[i]) * theta

    p[i]     <- 0.001 + 0.998 * p_raw[i]
    p_raw[i] <- alpha[i] * (1 - exp(-(gamma[i] * x[i])))

    logit(alpha[i]) <- par_alpha0[inds[i]] +
                        rw_alpha_pop[cohort_idx[i]] +
                        study_re[study_id[i]]

    log(gamma[i])   <- par_gamma0[inds[i]] +
                        rw_gamma_pop[cohort_idx[i]] 
  }

  log_theta ~ dnorm(3, 1)
  theta     <- exp(log_theta)

  ## shared/population cohort trajectories, i.e. trends over time in both prevalence & rates of infection by age
  rw_alpha_pop[1] ~ dnorm(0, rw_init_tau)
  rw_gamma_pop[1] ~ dnorm(0, rw_init_tau)
  for (t in 2:n_cohorts){
    rw_alpha_pop[t] ~ dnorm(rw_alpha_pop[t-1], tau_rw_alpha_pop)
    rw_gamma_pop[t] ~ dnorm(rw_gamma_pop[t-1], tau_rw_gamma_pop)
  }
  tau_rw_alpha_pop ~ dgamma(2, 0.5)
  tau_rw_gamma_pop ~ dgamma(2, 0.5)
  sd_rw_alpha_pop  <- 1/sqrt(tau_rw_alpha_pop)
  sd_rw_gamma_pop  <- 1/sqrt(tau_rw_gamma_pop)


  ## study-level random effects for alpha
  for (s in 1:n_studies){
    study_re[s] ~ dnorm(0, tau_study)
  }
  tau_study ~ dgamma(2, 0.5)
  sd_study  <- 1/sqrt(tau_study)

  ## race-specific intercepts (baseline prevalence and rates of infection at the referent cohort)
  for (j in 1:n_races){
    par_alpha0[j] ~ dt(mub[1], taub_diag[1], 1)
    par_gamma0[j] ~ dt(mub[3], taub_diag[3], 3)
  }
  mub[1] ~ dt(0, 1, 1)
  mub[3] ~ dt(-3.9, 1, 3)
  taub_diag[1] ~ dgamma(2, 1)
  taub_diag[3] ~ dgamma(2, 1)

  ForMCMCConv[1] <- mub[1]
  ForMCMCConv[2] <- mub[3]
  ForMCMCConv[3] <- theta
  ForMCMCConv[4] <- sd_rw_alpha_pop
  ForMCMCConv[5] <- sd_rw_gamma_pop
  ForMCMCConv[6] <- sd_study
}")


model_RW_expanded_continuous_bc <- paste0("model{

  ## Beta-binomial likelihood
  for (i in 1:N){
    y[i]   ~ dbinom(pi[i], n[i])
    pi[i]  ~ dbeta(a[i], b[i])
    a[i]   <- p[i] * theta
    b[i]   <- (1 - p[i]) * theta

    p[i]     <- 0.001 + 0.998 * p_raw[i]
    p_raw[i] <- alpha[i] * (1 - exp(-(gamma[i] * x[i])))

    logit(alpha[i]) <- par_alpha0[inds[i]] +
                        rw_alpha_pop[cohort_idx[i]] +
                        rw_alpha_dev[inds[i], cohort_idx[i]] +
                        study_re[study_id[i]]

    log(gamma[i])   <- par_gamma0[inds[i]] +
                        rw_gamma_pop[cohort_idx[i]] +
                        rw_gamma_dev[inds[i], cohort_idx[i]]
  }

  log_theta ~ dnorm(3, 1)
  theta     <- exp(log_theta)

  rw_alpha_pop[1] ~ dnorm(0, rw_init_tau)
  rw_gamma_pop[1] ~ dnorm(0, rw_init_tau)
  for (t in 2:n_cohorts){
    rw_alpha_pop[t] ~ dnorm(rw_alpha_pop[t-1], tau_rw_alpha_pop)
    rw_gamma_pop[t] ~ dnorm(rw_gamma_pop[t-1], tau_rw_gamma_pop)
  }
  tau_rw_alpha_pop ~ dgamma(2, 0.5)
  tau_rw_gamma_pop ~ dgamma(2, 0.5)
  sd_rw_alpha_pop  <- 1/sqrt(tau_rw_alpha_pop)
  sd_rw_gamma_pop  <- 1/sqrt(tau_rw_gamma_pop)

  ### race-specific deviations from the shared trajectory 
  for (j in 1:n_races){
   rw_alpha_dev[j,1] ~ dnorm(0, tau_rw_dev_alpha)
   rw_gamma_dev[j,1] ~ dnorm(0, tau_rw_dev_gamma)
   for (t in 2:n_cohorts){
    rw_alpha_dev[j,t] ~ dnorm(rw_alpha_dev[j,t-1], tau_rw_dev_alpha)
    rw_gamma_dev[j,t] ~ dnorm(rw_gamma_dev[j,t-1], tau_rw_dev_gamma)
   }
  }
  tau_rw_dev_alpha ~ dgamma(2, 1)
  tau_rw_dev_gamma ~ dgamma(2, 1)
  sd_rw_dev_alpha  <- 1/sqrt(tau_rw_dev_alpha)
  sd_rw_dev_gamma  <- 1/sqrt(tau_rw_dev_gamma)

  for (s in 1:n_studies){
    study_re[s] ~ dnorm(0, tau_study)
  }
  tau_study ~ dgamma(2, 0.5)
  sd_study  <- 1/sqrt(tau_study)

  for (j in 1:n_races){
    par_alpha0[j] ~ dt(mub[1], taub_diag[1], 1)
    par_gamma0[j] ~ dt(mub[3], taub_diag[3], 3)
  }
  mub[1] ~ dt(0, 1, 1)
  mub[3] ~ dt(-3.9, 1, 3)
  taub_diag[1] ~ dgamma(2, 1)
  taub_diag[3] ~ dgamma(2, 1)

  ForMCMCConv[1] <- mub[1]
  ForMCMCConv[2] <- mub[3]
  ForMCMCConv[3] <- theta
  ForMCMCConv[4] <- sd_rw_alpha_pop
  ForMCMCConv[5] <- sd_rw_gamma_pop
  ForMCMCConv[6] <- sd_rw_dev_alpha
  ForMCMCConv[7] <- sd_rw_dev_gamma
  ForMCMCConv[8] <- sd_study
}")

# *****************************************************************************
# Registry of all models
# *****************************************************************************
MODEL_LIST <- list(

  base_continuous_bc = list(
    model_string = model_base_continuous_bc,
    data         = base_data_continuous,
    monitor      =  c("ForMCMCConv",
               "mub", "var_b",
               "par", "study_re", "sd_study", "theta",
               "deviance")
  ),

  RW_continuous_bc = list(
    model_string = model_RW_continuous_bc,
    data         = rw_data,
    monitor      = c("ForMCMCConv", "par_alpha0", "par_gamma0",
               "rw_alpha_pop", "rw_gamma_pop",
               "study_re", "sd_rw_alpha_pop", "sd_rw_gamma_pop",
                "sd_study",
               "theta", "deviance")
  ),

  RW_expanded_continuous_bc = list(
    model_string = model_RW_expanded_continuous_bc,
    data         = rw_data,
    monitor      = c("ForMCMCConv", "par_alpha0", "par_gamma0",
               "rw_alpha_pop", "rw_gamma_pop", "rw_alpha_dev", "rw_gamma_dev",
               "study_re", "sd_rw_alpha_pop", "sd_rw_gamma_pop",
               "sd_rw_dev_alpha", "sd_rw_dev_gamma",
                "sd_study",
               "theta", "deviance")
  )

  
)

# *****************************************************************************
#### MCMC settings #### 
# *****************************************************************************
n_chains <- 3
n_adapt  <- 10000
n_burn   <- 500000
n_post   <- 20000
thinning <- 20

# per-chain RNG seeds AND over-dispersed starting values.
inits_list <- list(
  list(.RNG.name = "base::Wichmann-Hill",        .RNG.seed = 1, log_theta = 2.5),
  list(.RNG.name = "base::Marsaglia-Multicarry", .RNG.seed = 2, log_theta = 3.0),
  list(.RNG.name = "base::Super-Duper",          .RNG.seed = 3, log_theta = 3.5)
)


# *****************************************************************************
#### Function to run and extract all results & convergence metrics by model #### 
# *****************************************************************************
fit_one_model <- function(model_name, spec) {

  cat(sprintf("\n===== Fitting model: %s =====\n", model_name))
  t0 <- Sys.time()

  fit <- run.jags(
    model      = spec$model_string,
    data       = spec$data,
    inits      = inits_list,
    monitor    = union(spec$monitor, "deviance"),
    n.chains   = n_chains,
    adapt      = n_adapt,
    burnin     = n_burn,
    sample     = n_post,
    thin       = thinning,
    method     = "parallel",
    summarise  = FALSE
  )

  runtime_mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("Finished %s in %.1f minutes (adapt=%d, burnin=%d, sample=%d, thin=%d).\n",
              model_name, runtime_mins, n_adapt, n_burn, n_post, thinning))

  mcmc_all <- as.mcmc.list(fit)

  # --- Convergence: ForMCMCConv block ---
  conv_cols <- grep("^ForMCMCConv", varnames(mcmc_all))
  mcmc_conv <- mcmc_all[, conv_cols, drop = FALSE]
  gelman_conv <- tryCatch(gelman.diag(mcmc_conv, multivariate = FALSE),
                           error = function(e) NULL)

  # --- Convergence: par[] block (race-level REs) ---
  par_cols <- grep("^par\\[|^par_alpha0|^par_gamma0", varnames(mcmc_all))
  mcmc_par <- mcmc_all[, par_cols, drop = FALSE]
  gelman_par <- tryCatch(gelman.diag(mcmc_par, multivariate = FALSE),
                          error = function(e) NULL)

  mean_gelman_conv <- if (!is.null(gelman_conv)) mean(gelman_conv$psrf[, 1]) else NA
  mean_gelman_par  <- if (!is.null(gelman_par))  mean(gelman_par$psrf[, 1])  else NA
  mean_gelman_all  <- mean(c(mean_gelman_conv, mean_gelman_par), na.rm = TRUE)

  max_gelman_conv <- if (!is.null(gelman_conv)) max(gelman_conv$psrf[, 1]) else NA
  max_gelman_par  <- if (!is.null(gelman_par))  max(gelman_par$psrf[, 1])  else NA
  max_gelman_all  <- suppressWarnings(max(c(max_gelman_conv, max_gelman_par), na.rm = TRUE))
  if (!is.finite(max_gelman_all)) max_gelman_all <- NA_real_

  n_above_1.1 <- sum(c(if (!is.null(gelman_conv)) gelman_conv$psrf[, 1],
                       if (!is.null(gelman_par))  gelman_par$psrf[, 1]) > 1.1, na.rm = TRUE)

  # flag as not converged if max PSRF > 1.1 
  if (!is.na(max_gelman_all) && max_gelman_all > 1.1) {
    warning("Model '", model_name, "': maximum PSRF is ", round(max_gelman_all, 3),
            " over ", n_above_1.1, " parameter(s) above 1.1. Chains have NOT ",
            "converged; do not report these results.",
            call. = FALSE, immediate. = TRUE)
  }

  # --- DIC via variance-of-deviance pD approximation ---
  dev_col <- which(varnames(mcmc_all) == "deviance")
  deviance_draws <- as.matrix(mcmc_all)[, dev_col]
  deviance_mean  <- mean(deviance_draws)
  pD_mean        <- var(deviance_draws) / 2
  dic_val        <- deviance_mean + pD_mean

  result <- list(
    model_name       = model_name,
    fit              = fit,
    mcmc_all         = mcmc_all,
    gelman_conv      = gelman_conv,
    gelman_par       = gelman_par,
    mean_gelman_conv = mean_gelman_conv,
    mean_gelman_par  = mean_gelman_par,
    mean_gelman_all  = mean_gelman_all,
    max_gelman_all   = max_gelman_all,
    n_PSRF_above_1.1 = n_above_1.1,
    deviance         = deviance_mean,
    pD               = pD_mean,
    dic              = dic_val,
    runtime_mins     = runtime_mins,
    n_chains         = n_chains,
    n_adapt          = n_adapt,
    n_burn           = n_burn,
    n_post           = n_post,
    thin             = thinning
  )

  saveRDS(result, file.path(out_dir, paste0("fit_", model_name, ".rds")))

  summary_row <- data.frame(
    model            = model_name,
    DIC              = dic_val,
    deviance         = deviance_mean,
    pD               = pD_mean,
    mean_PSRF_conv   = mean_gelman_conv,
    mean_PSRF_par    = mean_gelman_par,
    mean_PSRF_all    = mean_gelman_all,
    max_PSRF_all     = max_gelman_all,
    n_PSRF_above_1.1 = n_above_1.1,
    converged        = !is.na(max_gelman_all) & max_gelman_all <= 1.1,
    runtime_mins     = runtime_mins,
    n_burn           = n_burn,
    n_post           = n_post,
    thin             = thinning,
    timestamp        = as.character(Sys.time())
  )

  write_csv(summary_row,
            file.path(out_dir, paste0("summary_", model_name, ".csv")))

  cat(sprintf("Model: %-26s DIC = %.1f | mean PSRF = %.4f | max PSRF = %.4f%s\n",
              model_name, dic_val, mean_gelman_all, max_gelman_all,
              if (!is.na(max_gelman_all) && max_gelman_all > 1.1) "  *** NOT CONVERGED ***" else ""))

  summary_row
}

# *****************************************************************************
#### Run one model or default of all models ####
# *****************************************************************************
args <- commandArgs(trailingOnly = TRUE)
model_names <- names(MODEL_LIST)

if (length(args) >= 1) {
  run_names <- intersect(args, model_names)
  if (length(run_names) == 0) {
    stop("Unrecognized model name(s). Available: ", paste(model_names, collapse = ", "))
  }
} else {
  run_names <- model_names
}

cat("Running model(s): ", paste(run_names, collapse = ", "), "\n")

all_summaries <- purrr::map_dfr(run_names, function(nm) {
  fit_one_model(nm, MODEL_LIST[[nm]])
})

# *****************************************************************************
#### Combined comparison table ####
# *****************************************************************************
existing_csvs <- list.files(out_dir, pattern = "^summary_.*\\.csv$", full.names = TRUE)
combined <- purrr::map_dfr(existing_csvs, read_csv, show_col_types = FALSE) %>%
  distinct(model, .keep_all = TRUE) %>%
  arrange(DIC)

write_csv(combined, file.path(out_dir, "model_comparison_sensitivity_analysis_ALL.csv"))

cat("\n===== Model comparison (sorted by DIC) =====\n")
print(combined)

cat(sprintf("\nAll results saved to: %s\n", normalizePath(out_dir)))
