
# *****************************************************************************
#### Load necessary packages ####
# *****************************************************************************

library(tidyverse)
library(data.table)
library(coda)
library(matrixStats)
library(gt)

dir.create("figures_sep22_2026", showWarnings = FALSE)
dir.create("figures_sep22_2026/continuous_cohort_figs", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE)

# *****************************************************************************
#### Shared data prep (identical across all three models) ####
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

ilogit <- function(logodds) 1 / (1 + exp(-logodds))

n_races     <- 5
fin_age     <- 90
ages_vec    <- 0:fin_age
n_age       <- fin_age + 1
race_index  <- 1:5
race_labels <- c("Hispanic", "AI/AN", "API", "NHB", "NHW")
race_short_labels <- c("H", "AIAN", "NHA", "NHB", "NHW")

ref_ages <- c(5, 10, fin_age)

# RW cohort indexing for MODEL 2 and MODEL 3:
#   "observed" = one RW node per observed birth year (fit_RW_*continuous_bc)
#   "fullyear" = one RW node per calendar year, min to max observed cohort
#                (fit_RW_*fullyear_continuous_bc); years with no data are latent
# Sets which fits are read, the cohort lookup, and a suffix on every RW output file.
RW_INDEXING <- "fullyear"
stopifnot(RW_INDEXING %in% c("observed", "fullyear"))
rw_fit_tag <- if (RW_INDEXING == "fullyear") "fullyear_" else ""
rw_out_tag <- if (RW_INDEXING == "fullyear") "_fullyear" else ""

# *****************************************************************************
#### Shared helper functions ####
# *****************************************************************************

# maps RW node index -> birth cohort year, matching the indexing in
# all_models_sensitivity_analysis_cluster.R; observed = year has data
make_cohort_lookup <- function() {
  obs_years <- sort(unique(round(data_model_filtered$birth_cohort)))
  years <- if (RW_INDEXING == "fullyear") seq(min(obs_years), max(obs_years)) else obs_years
  data.frame(cohort_idx = seq_along(years), birth_cohort = years,
             observed = years %in% obs_years)
}

get_vec <- function(mat_all, base_name, idx) {
  mat_all[, paste0(base_name, "[", idx, "]")]
}

# Effective sample size, adjusted for beta-binomial overdispersion (theta), for
# the selected model only (RW_continuous_bc -- Appendix Section D).
#
# For y ~ BetaBinomial(n, p*theta, (1-p)*theta), Var(y) is inflated relative to
# a plain Binomial(n, p) by a design effect DEFF = (n + theta) / (1 + theta).
# The effective sample size is the plain-binomial n that would give the same
# precision: n_eff = n / DEFF = n*(1 + theta) / (n + theta), computed per row
# and per posterior draw of theta, then summed across rows -- same structure
# as total_pop, but downweighted for between-study/subgroup overdispersion.
compute_effective_n <- function(fit_path, n_vec) {
  fit_obj <- readRDS(fit_path)
  theta_draws <- as.matrix(fit_obj$mcmc_all)[, "theta"]
  eff_n_by_draw <- colSums(outer(n_vec, theta_draws, function(n, th) n * (1 + th) / (n + th)))
  data.frame(
    total_pop  = sum(n_vec),
    eff_n_mean = mean(eff_n_by_draw),
    eff_n_LB   = quantile(eff_n_by_draw, 0.025),
    eff_n_UB   = quantile(eff_n_by_draw, 0.975)
  )
}

to_long_stat <- function(mat, stat_name, idx_name) {
  dt <- as.data.table(reshape2::melt(mat, varnames = c("age_idx", idx_name),
                                      value.name = stat_name))
  dt
}

merge_stats_long <- function(mean_mat, sd_mat, lb_mat, ub_mat, idx_name) {
  dt <- to_long_stat(mean_mat, "Mean", idx_name)[
        to_long_stat(sd_mat,   "SD",   idx_name), on = c("age_idx", idx_name)][
        to_long_stat(lb_mat,   "LB",   idx_name), on = c("age_idx", idx_name)][
        to_long_stat(ub_mat,   "UB",   idx_name), on = c("age_idx", idx_name)]
  dt
}

run_predictions <- function(par_fn, grid_df, idx_name, race_index, race_labels,
                             ages_vec = 0:90) {
  n_age  <- length(ages_vec)
  n_grid <- nrow(grid_df)

  prev_list <- vector("list", length(race_index))
  foi_list  <- vector("list", length(race_index))

  for (r in race_index) {

    prev_mean <- prev_sd <- prev_lb <- prev_ub <- matrix(NA_real_, n_age, n_grid)
    foi_mean  <- foi_sd  <- foi_lb  <- foi_ub  <- matrix(NA_real_, n_age, n_grid)

    for (t in seq_len(n_grid)) {

      ag <- par_fn(r, t)
      a_vec <- ag$a_vec
      g_vec <- ag$g_vec

      g_age         <- outer(ages_vec, g_vec, FUN = "*")
      exp_neg       <- exp(-g_age)
      one_minus_exp <- 1 - exp_neg

      a_rep <- rep(a_vec, each = n_age)
      g_rep <- rep(g_vec, each = n_age)

      p_raw_mat <- one_minus_exp * a_rep
      p_mat     <- 0.001 + 0.998 * p_raw_mat

      # force of infection: lambda(a) = (dF/da) / (1 - F(a))   [Appendix Eq. 3]
      foi_num <- a_rep * g_rep * exp_neg
      foi_den <- pmax(1 - p_raw_mat, 1e-10)
      foi_mat <- foi_num / foi_den
      foi_mat[ages_vec == 0, ] <- NA

      prev_mean[, t] <- rowMeans2(p_raw_mat, na.rm = TRUE)
      prev_sd[, t]   <- rowSds(p_raw_mat,   na.rm = TRUE)
      pq <- rowQuantiles(p_raw_mat, probs = c(0.025, 0.975), na.rm = TRUE)
      prev_lb[, t] <- pq[, 1]; prev_ub[, t] <- pq[, 2]

      foi_mean[, t] <- rowMeans2(foi_mat, na.rm = TRUE)
      foi_sd[, t]   <- rowSds(foi_mat,   na.rm = TRUE)
      fq <- rowQuantiles(foi_mat, probs = c(0.025, 0.975), na.rm = TRUE)
      foi_lb[, t] <- fq[, 1]; foi_ub[, t] <- fq[, 2]
    }

    dt_prev <- merge_stats_long(prev_mean, prev_sd, prev_lb, prev_ub, idx_name)
    dt_prev[, `:=`(age = age_idx - 1, race = race_labels[r], type = "Prevalence")]
    dt_prev <- merge(dt_prev, grid_df, by = idx_name, all.x = TRUE)

    dt_foi <- merge_stats_long(foi_mean, foi_sd, foi_lb, foi_ub, idx_name)
    dt_foi[, `:=`(age = age_idx - 1, race = race_labels[r], type = "Force of Infection")]
    dt_foi <- merge(dt_foi, grid_df, by = idx_name, all.x = TRUE)

    prev_list[[r]] <- dt_prev
    foi_list[[r]]  <- dt_foi

    cat("Finished predictions for race:", race_labels[r], "\n")
  }

  list(prev = bind_rows(prev_list), foi = bind_rows(foi_list))
}

run_aoi <- function(gamma_fn, grid_df, idx_name, race_index, race_labels, n_iter) {
  n_grid <- nrow(grid_df)
  aoi_list <- vector("list", length(race_index))

  for (r in race_index) {
    aoi_mat <- matrix(NA, nrow = n_iter, ncol = n_grid)

    for (t in seq_len(n_grid)) {
      for (iter in 1:n_iter) {
        g <- gamma_fn(r, t, iter)
        aoi_mat[iter, t] <- log(2) / g
      }
    }

    dt_aoi <- data.frame(
      idx  = seq_len(n_grid),
      Mean = apply(aoi_mat, 2, mean),
      LB   = apply(aoi_mat, 2, quantile, 0.025),
      UB   = apply(aoi_mat, 2, quantile, 0.975)
    )
    names(dt_aoi)[1] <- idx_name
    dt_aoi <- dt_aoi %>%
      left_join(grid_df, by = idx_name) %>%
      mutate(race = race_labels[r], type = "Median age of Infection")

    aoi_list[[r]] <- dt_aoi
    cat("Finished median age of infection for race:", race_labels[r], "\n")
  }

  bind_rows(aoi_list)
}

generate_traceplot_function <- function(draws, labels, ncol_facet = 1, n_chains = 3) {
  draws_df <- as.data.frame(draws)
  colnames(draws_df) <- labels
  n_iter_per_chain <- nrow(draws_df) / n_chains
  stopifnot(n_iter_per_chain == round(n_iter_per_chain))
  draws_df$iteration <- rep(seq_len(n_iter_per_chain), n_chains)
  draws_df$chain <- rep(factor(seq_len(n_chains)), each = n_iter_per_chain)

  draws_df_long <- pivot_longer(draws_df,
                                 cols = -c(iteration,chain),
                                 names_to = "parameter",
                                 values_to = "value")

  ggplot(draws_df_long, aes(x = iteration, y = value, color = chain)) +
    geom_line(alpha = 0.7,size=1) +
    facet_wrap(~ parameter, scales = "free_y", ncol = ncol_facet) +
    theme_minimal() +
    labs(x = "Iteration", y = "Estimate") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
          strip.text = element_text(size = 20),
        legend.text=element_text(size=20), legend.title=element_text(size=20))
}

generate_densplot_function <- function(draws, labels, ncol_facet = 1) {
  draws_df <- as.data.frame(draws)
  colnames(draws_df) <- labels
  draws_df$iteration <- seq_len(nrow(draws_df))

  draws_df_long <- pivot_longer(draws_df,
                                 cols = -iteration,
                                 names_to = "parameter",
                                 values_to = "value")

  ggplot(draws_df_long, aes(x = value, fill = parameter)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ parameter, scales = "free_y", ncol = ncol_facet) +
    theme_bw() +
    labs(x = "Estimate", y = "Density") +
    theme(legend.position = "none",
           axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
          strip.text = element_text(size = 20))
}

make_max_prev_by_race_cohort <- function(df_p_hat_race_final, ref_age, filter_years = NULL) {
  out <- df_p_hat_race_final %>% filter(age == ref_age)
  if (!is.null(filter_years)) out <- out %>% filter(birth_cohort %in% filter_years)
  out %>% mutate(ref_age = ref_age)
}

# *****************************************************************************************************
# *****************************************************************************************************
#### MODEL 1: Base continuous (bc_scaled as a continuous linear covariate, no random walk) ####
# *****************************************************************************************************
# *****************************************************************************************************

cluster_fit <- readRDS("./results_cluster_sensitivity_analysis/fit_base_continuous_bc.rds")
mcmc_all    <- cluster_fit$mcmc_all
mat_all     <- as.matrix(mcmc_all)

n_iter <- nrow(mat_all)

bc_range    <- range(data_model_filtered$bc_scaled, na.rm = TRUE)
n_bc_points <- 67

bc_scaled_grid     <- seq(bc_range[1], bc_range[2], length.out = n_bc_points)
birth_cohort_grid  <- bc_scaled_grid + min(data_model_filtered$birth_cohort)

bc_lookup <- data.frame(
  bc_point_idx = seq_len(n_bc_points),
  bc_scaled    = bc_scaled_grid,
  birth_cohort = round(birth_cohort_grid)
)

get_par_matrix <- function(mat_all, param_col, n_races) {
  cols <- paste0("par[", 1:n_races, ",", param_col, "]")
  mat_all[, cols, drop = FALSE]
}

m_alpha0 <- t(get_par_matrix(mat_all, 1, n_races))   # [race x iter]
m_alpha1 <- t(get_par_matrix(mat_all, 2, n_races))
m_gamma0 <- t(get_par_matrix(mat_all, 3, n_races))
m_gamma1 <- t(get_par_matrix(mat_all, 4, n_races))

n_studies <- length(unique(as.numeric(as.factor(data_model_filtered$study_id))))
study_re_draws <- t(sapply(1:n_studies, function(s) get_vec(mat_all, "study_re", s)))
theta_draws    <- mat_all[, "theta"]
sd_study_draws <- mat_all[, "sd_study"]

par_fn_base <- function(r, t) {
  bc_val <- bc_scaled_grid[t]
  list(
    a_vec = ilogit(m_alpha0[r, ] + m_alpha1[r, ] * bc_val),
    g_vec = exp(   m_gamma0[r, ] + m_gamma1[r, ] * bc_val)
  )
}

pred_base <- run_predictions(par_fn_base, bc_lookup, "bc_point_idx", race_index, race_labels, ages_vec)
df_p_hat_race_final  <- pred_base$prev
df_hp_hat_race_final <- pred_base$foi

for (r in race_index) {
  dt_prev <- df_p_hat_race_final  %>% filter(race == race_labels[r])
  dt_foi  <- df_hp_hat_race_final %>% filter(race == race_labels[r])
  save(dt_prev, dt_foi,
       file = paste0("results/base_continuous_dt_prev_hat_race_", r, "_",
                      race_short_labels[r], "_BaseContinuous.Rdata"))
}

df_foi_curves_subset <- df_hp_hat_race_final %>% filter(age > 0, Mean < 1)

ggplot(df_foi_curves_subset,
       aes(x = age, y = Mean, color = birth_cohort, group = birth_cohort)) +
  geom_line(linewidth = 1) +
  facet_wrap(~race, ncol = 3) +
  scale_color_viridis_c(name = "Birth cohort\nyear") +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  xlab("Age") + ylab("Predicted force of infection") +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "bold",size=20))

ggsave("./figures_sep22_2026/continuous_cohort_figs/predicted_foi_by_cohort_year_BaseContinuous.png", width = 12, height = 8)

df_p_curves_subset <- df_p_hat_race_final %>% filter(age > 0)

ggplot(df_p_curves_subset,
       aes(x = age, y = Mean, color = birth_cohort, group = birth_cohort)) +
  geom_line(linewidth = 1) +
  facet_wrap(~race, ncol = 3) +
  scale_color_viridis_c(name = "Birth cohort\nyear") +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  xlab("Age") + ylab("Predicted seroprevalence") +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "bold",size=20))

ggsave("./figures_sep22_2026/continuous_cohort_figs/predicted_prev_by_cohort_year_BaseContinuous.png", width = 12, height = 8)

plot_max_prev_base <- function(df, ref_age) {
  df <- df%>%
    mutate(Mean=Mean*100)%>%
    mutate(LB=LB*100)%>%
    mutate(UB=UB*100)
  ggplot(df,
         aes(x = race, y = Mean, color = factor(birth_cohort),
             ymin = LB, ymax = UB, group = factor(birth_cohort))) +
    geom_pointrange(position = position_dodge(width = 0.6), size = 0.7) +
    scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 10)) +
    labs(x = "", y = paste0("Seroprevalence (age ", ref_age, ")"), color = "Birth cohort\nyear") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}

max_prev_by_race_cohort_list <- lapply(ref_ages, function(a)
  make_max_prev_by_race_cohort(df_p_hat_race_final, a))
names(max_prev_by_race_cohort_list) <- paste0("age_", ref_ages)
max_prev_by_race_cohort <- max_prev_by_race_cohort_list[[paste0("age_", fin_age)]]

for (ref_age in ref_ages) {
  df_this <- max_prev_by_race_cohort_list[[paste0("age_", ref_age)]]
  p_this  <- plot_max_prev_base(df_this, ref_age)

  out_file <- if (ref_age == fin_age) {
    "./figures_sep22_2026/continuous_cohort_figs/max_prev_by_race_cohort_BaseContinuous.png"
  } else {
    paste0("./figures_sep22_2026/continuous_cohort_figs/max_prev_by_race_cohort_age", ref_age, "_BaseContinuous.png")
  }
  ggsave(out_file, p_this, width = 11, height = 8)
}

max_prev_by_race_cohort_all_ages <- bind_rows(max_prev_by_race_cohort_list)
save(max_prev_by_race_cohort_all_ages,
     file = "results/base_continuous_max_prev_age5_10_90.Rdata")

## median age of infection
gamma_fn_base <- function(r, t, iter) {
  bc_val <- bc_scaled_grid[t]
  exp(m_gamma0[r, iter] + m_gamma1[r, iter] * bc_val)
}
df_aoi_final <- run_aoi(gamma_fn_base, bc_lookup, "bc_point_idx", race_index, race_labels, n_iter)

for (r in race_index) {
  dt_aoi <- df_aoi_final %>% filter(race == race_labels[r])
  save(dt_aoi, file = paste0("results/base_continuous_dt_aoi_race_", r, "_",
                              race_short_labels[r], "_BaseContinuous.Rdata"))
}

write_csv(df_aoi_final, "results/median_age_inf_BaseContinuous.csv")

## convergence plots
alpha0_labels <- paste0("alpha intercept ", race_labels)
alpha1_labels <- paste0("alpha slope ", race_labels)
gamma0_labels <- paste0("gamma intercept ", race_labels)
gamma1_labels <- paste0("gamma slope ", race_labels)

generate_traceplot_function(t(m_alpha0), alpha0_labels, ncol_facet = 1)
ggsave("./figures_sep22_2026/continuous_cohort_figs/alpha0_traceplot_BaseContinuous.png", width = 12, height = 8)
generate_traceplot_function(t(m_alpha1), alpha1_labels, ncol_facet = 1)
ggsave("./figures_sep22_2026/continuous_cohort_figs/alpha1_traceplot_BaseContinuous.png", width = 12, height = 8)
generate_traceplot_function(t(m_gamma0), gamma0_labels, ncol_facet = 1)
ggsave("./figures_sep22_2026/continuous_cohort_figs/gamma0_traceplot_BaseContinuous.png", width = 12, height = 8)
generate_traceplot_function(t(m_gamma1), gamma1_labels, ncol_facet = 1)
ggsave("./figures_sep22_2026/continuous_cohort_figs/gamma1_traceplot_BaseContinuous.png", width = 12, height = 8)

generate_densplot_function(t(m_alpha0), alpha0_labels, ncol_facet = 1)
ggsave("./figures_sep22_2026/continuous_cohort_figs/alpha0_densplot_BaseContinuous.png", width = 12, height = 8)
generate_densplot_function(t(m_alpha1), alpha1_labels, ncol_facet = 1)
ggsave("./figures_sep22_2026/continuous_cohort_figs/alpha1_densplot_BaseContinuous.png", width = 12, height = 8)
generate_densplot_function(t(m_gamma0), gamma0_labels, ncol_facet = 1)
ggsave("./figures_sep22_2026/continuous_cohort_figs/gamma0_densplot_BaseContinuous.png", width = 12, height = 8)
generate_densplot_function(t(m_gamma1), gamma1_labels, ncol_facet = 1)
ggsave("./figures_sep22_2026/continuous_cohort_figs/gamma1_densplot_BaseContinuous.png", width = 12, height = 8)

gelman_conv <- cluster_fit$gelman_conv
gelman_multivariate_PSRF <- as.data.frame(gelman_conv$psrf[, 1])
colnames(gelman_multivariate_PSRF) <- "PSRF"
gelman_multivariate_PSRF$param <- c(
  "mu_alpha0", "mu_alpha1", "mu_gamma0", "mu_gamma1",
  "var_alpha0", "var_alpha1", "var_gamma0", "var_gamma1",
  "cov_alpha0_1", "cov_alpha0_gamma0",
  "cov_gamma0_alpha1", "cov_gamma1_alpha0",
  "cov_gamma1_alpha1", "cov_gamma0_1",
  "overdispersion_parameter", "sd_study"
)
# gtsave() to PNG requires the webshot2 + chromote packages and a headless
# Chrome install. When those are missing it throws, which previously aborted
# the whole per-model loop after the first model -- losing every subsequent
# model's figures for the sake of a diagnostics table. This always writes the
# table as CSV (which is what you actually want for reporting PSRF values) and
# treats the PNG as best-effort.
save_psrf_table <- function(df, png_path) {
  csv_path <- sub("\\.png$", ".csv", png_path)
  readr::write_csv(df, csv_path)
  ok <- requireNamespace("webshot2", quietly = TRUE) &&
        requireNamespace("chromote", quietly = TRUE)
  if (!ok) {
    message("  note: wrote ", csv_path,
            " (PNG skipped; install webshot2 + chromote for image output)")
    return(invisible(FALSE))
  }
  res <- try(gt::gtsave(gt::gt(df), png_path), silent = TRUE)
  if (inherits(res, "try-error")) {
    message("  note: wrote ", csv_path, " but gtsave() to PNG failed: ",
            conditionMessage(attr(res, "condition")))
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

save_psrf_table(cbind.data.frame(parameter = gelman_multivariate_PSRF$param,
                                 PSRF = gelman_multivariate_PSRF$PSRF),
                "./figures_sep22_2026/continuous_cohort_figs/multivariate_PSRF_BaseContinuous.png")

gelman_par <- cluster_fit$gelman_par
gelman_univariate_PSRF <- as.data.frame(gelman_par$psrf[, 1])
colnames(gelman_univariate_PSRF) <- "PSRF"
gelman_univariate_PSRF$param <- c(
  "log_odds_alpha_H", "log_odds_alpha_AIAN", "log_odds_alpha_API",
  "log_odds_alpha_NHB", "log_odds_alpha_NHW",
  "delta_log_odds_alpha_H", "delta_log_odds_alpha_AIAN", "delta_log_odds_alpha_API",
  "delta_log_odds_alpha_NHB", "delta_log_odds_alpha_NHW",
  "log_gamma_H", "log_gamma_AIAN", "log_gamma_API",
  "log_gamma_NHB", "log_gamma_NHW",
  "delta_log_gamma_H", "delta_log_gamma_AIAN", "delta_log_gamma_API",
  "delta_log_gamma_NHB", "delta_log_gamma_NHW"
)
save_psrf_table(cbind.data.frame(parameter = gelman_univariate_PSRF$param,
                                 PSRF = gelman_univariate_PSRF$PSRF),
                "./figures_sep22_2026/continuous_cohort_figs/univariate_PSRF_BaseContinuous.png")


# *****************************************************************************************************
# *****************************************************************************************************
#### MODEL 2: RW cohort (continuous birth cohort, random-walk cohort effect) ####
# *****************************************************************************************************
# *****************************************************************************************************

cluster_fit <- readRDS(paste0("./results_cluster_sensitivity_analysis/fit_RW_", rw_fit_tag, "continuous_bc.rds"))
mcmc_all    <- cluster_fit$mcmc_all
mat_all     <- as.matrix(mcmc_all)

n_iter     <- nrow(mat_all)
cohort_lookup <- make_cohort_lookup()
n_cohorts     <- nrow(cohort_lookup)

par_alpha0_draws <- t(sapply(1:n_races, function(j) get_vec(mat_all, "par_alpha0", j)))
par_gamma0_draws <- t(sapply(1:n_races, function(j) get_vec(mat_all, "par_gamma0", j)))

rw_alpha_pop_draws <- t(sapply(1:n_cohorts, function(t) get_vec(mat_all, "rw_alpha_pop", t)))
rw_gamma_pop_draws <- t(sapply(1:n_cohorts, function(t) get_vec(mat_all, "rw_gamma_pop", t)))

par_fn_rw <- function(r, t) {
  list(
    a_vec = ilogit(par_alpha0_draws[r, ] + rw_alpha_pop_draws[t, ]),
    g_vec = exp(   par_gamma0_draws[r, ] + rw_gamma_pop_draws[t, ])
  )
}

pred_rw <- run_predictions(par_fn_rw, cohort_lookup, "cohort_idx", race_index, race_labels, ages_vec)
df_p_hat_race_final  <- pred_rw$prev
df_hp_hat_race_final <- pred_rw$foi

for (r in race_index) {
  dt_prev <- df_p_hat_race_final  %>% filter(race == race_labels[r])
  dt_foi  <- df_hp_hat_race_final %>% filter(race == race_labels[r])
  save(dt_prev, dt_foi,
       file = paste0("results/rw_cohort_dt_prev_hat_race_", r, "_",
                      race_short_labels[r], paste0("_RWcohort", rw_out_tag, ".Rdata")))
}

## FOI plot
df_foi_curves_subset <- df_hp_hat_race_final %>% filter(age > 0, Mean < 1)
ggplot(df_foi_curves_subset,
       aes(x = age, y = Mean, color = birth_cohort, group = birth_cohort)) +
  geom_line(linewidth = 1) +
  facet_wrap(~race, ncol = 3) +
  scale_color_viridis_c(name = "Birth cohort\nyear") +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  xlab("Age") + ylab("Predicted force of infection") +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "bold", size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20, hjust = 1, vjust = 1, angle = 45),
        legend.title = element_text(size = 20), legend.position = "bottom")

ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/predicted_foi_by_cohort_year_RWcohort", rw_out_tag, ".png"),
       width = 12, height = 8)

## prevalence plot
df_p_curves_subset <- df_p_hat_race_final %>% filter(age > 0)
ggplot(df_p_curves_subset,
       aes(x = age, y = Mean, color = birth_cohort, group = birth_cohort)) +
  geom_line(linewidth = 1) +
  facet_wrap(~race, ncol = 3) +
  scale_color_viridis_c(name = "Birth cohort\nyear") +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  xlab("Age") + ylab("Predicted seroprevalence") +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "bold",size=20))

ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/predicted_prev_by_cohort_year_RWcohort", rw_out_tag, ".png"),
       width = 12, height = 8)

## max prevalence at multiple ages
max_prev_by_race_cohort_list <- lapply(ref_ages, function(a)
  make_max_prev_by_race_cohort(df_p_hat_race_final, a))
names(max_prev_by_race_cohort_list) <- paste0("age_", ref_ages)

plot_max_prev_rw <- function(df, ref_age) {
  df <- df%>%
    mutate(Mean=Mean*100)%>%
    mutate(LB=LB*100)%>%
    mutate(UB=UB*100)
   df <- df%>%
     mutate(race=ifelse(race=="NH American Indian","AI/AN",
    ifelse(race=="NH Asian","API",
  ifelse(race=="NH Black","NHB",
  ifelse(race=="NH White","NHW",
  "Hispanic")))))

  ggplot(df,
         aes(x = race, y = Mean, color = birth_cohort,
             ymin = LB, ymax = UB, group = birth_cohort)) +
    geom_pointrange(position = position_dodge(width = 0.6), size = 0.5) +
    scale_color_viridis_c(name = "Birth cohort\nyear") +
    scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
    labs(x = "", y = paste0("Seroprevalence (age ", ref_age, ") (%)"),
         color = "Birth cohort\nyear") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.position = "none",
          legend.text = element_text(size = 20, angle = 45, vjust = 1, hjust = 1),
          legend.title = element_text(size = 20))
}

for (ref_age in ref_ages) {
  df_this <- max_prev_by_race_cohort_list[[paste0("age_", ref_age)]]
  p_this  <- plot_max_prev_rw(df_this, ref_age)

  out_file <- if (ref_age == fin_age) {
    paste0("./figures_sep22_2026/continuous_cohort_figs/max_prev_by_race_cohort_RWcohort", rw_out_tag, ".png")
  } else {
    paste0("./figures_sep22_2026/continuous_cohort_figs/max_prev_by_race_cohort_age", ref_age, "_RWcohort", rw_out_tag, ".png")
  }
  ggsave(out_file, p_this, width = 11, height = 8)
}

max_prev_by_race_cohort_all_ages <- bind_rows(max_prev_by_race_cohort_list)
save(max_prev_by_race_cohort_all_ages,
     file = paste0("results/rw_cohort_max_prev_age5_10_90", rw_out_tag, ".Rdata"))

## median age of infection
gamma_fn_rw <- function(r, t, iter) exp(par_gamma0_draws[r, iter] + rw_gamma_pop_draws[t, iter])
df_aoi_rw <- run_aoi(gamma_fn_rw, cohort_lookup, "cohort_idx", race_index, race_labels, n_iter) %>%
  mutate(type = "Median age of Infection")

for (r in race_index) {
  dt_aoi <- df_aoi_rw %>% filter(race == race_labels[r])
  save(dt_aoi, file = paste0("results/rw_cohort_dt_aoi_race_", r, "_",
                              race_short_labels[r], paste0("_RWcohort", rw_out_tag, ".Rdata")))
}

write_csv(df_aoi_rw, paste0("results/median_age_inf_RW_cohort", rw_out_tag, ".csv"))

order_race_labels_aoi_rw=c("API", "NHW", "Hispanic", "NHB", "AI/AN") # ordered by magnitude of latest estimates, to correspond to binary cohort figures 
df_aoi_rw$race<-factor(df_aoi_rw$race, levels = order_race_labels_aoi_rw)

ggplot(df_aoi_rw,
       aes(x = race, y = Mean, #reorder(race, Mean, decreasing = TRUE)
           ymin = LB, ymax = UB,
           group = birth_cohort, color = birth_cohort)) +
  geom_point()+
  geom_pointrange(size = 0.3, position = position_dodge(width = 0.8)) +
  geom_errorbar(width = 0.3, linewidth = 0.5, position = position_dodge(width = 0.8)) +
  scale_color_viridis_c(name = "Birth cohort\nyear") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 10, angle = 45, hjust = 1, vjust = 1)) +
  ylab("Median age of infection (years)") +
  labs(color = "Birth cohort year")

ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/median_age_inf_RW_cohort", rw_out_tag, ".png"))

## convergence plots/tables
m_alpha0 <- mat_all[, paste0("par_alpha0[", 1:n_races, "]"), drop = FALSE]
m_gamma0 <- mat_all[, paste0("par_gamma0[", 1:n_races, "]"), drop = FALSE]
alpha0_labels <- paste0("alpha intercept ", race_labels)
gamma0_labels <- paste0("gamma intercept ", race_labels)

ref_cohort_bins  <- round(seq(1, n_cohorts, length.out = 6))
ref_cohort_years <- cohort_lookup$birth_cohort[ref_cohort_bins]

m_alpha_pop <- mat_all[, paste0("rw_alpha_pop[", ref_cohort_bins, "]"), drop = FALSE]
m_gamma_pop <- mat_all[, paste0("rw_gamma_pop[", ref_cohort_bins, "]"), drop = FALSE]
alpha_pop_labels <- paste0("alpha_pop, cohort=", ref_cohort_years)
gamma_pop_labels <- paste0("gamma_pop, cohort=", ref_cohort_years)

generate_traceplot_function(m_alpha0, alpha0_labels, ncol_facet = 1)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/alpha0_traceplot_RWcohort", rw_out_tag, ".png"), width = 12, height = 8)
generate_traceplot_function(m_gamma0, gamma0_labels, ncol_facet = 1)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/gamma0_traceplot_RWcohort", rw_out_tag, ".png"), width = 12, height = 8)
generate_traceplot_function(m_alpha_pop, alpha_pop_labels, ncol_facet = 2)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/alpha_pop_traceplot_RWcohort", rw_out_tag, ".png"), width = 14, height = 10)
generate_traceplot_function(m_gamma_pop, gamma_pop_labels, ncol_facet = 2)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/gamma_pop_traceplot_RWcohort", rw_out_tag, ".png"), width = 14, height = 10)

generate_densplot_function(m_alpha0, alpha0_labels, ncol_facet = 1)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/alpha0_densplot_RWcohort", rw_out_tag, ".png"), width = 12, height = 8)
generate_densplot_function(m_gamma0, gamma0_labels, ncol_facet = 1)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/gamma0_densplot_RWcohort", rw_out_tag, ".png"), width = 12, height = 8)
generate_densplot_function(m_alpha_pop, alpha_pop_labels, ncol_facet = 2)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/alpha_pop_densplot_RWcohort", rw_out_tag, ".png"), width = 14, height = 10)
generate_densplot_function(m_gamma_pop, gamma_pop_labels, ncol_facet = 2)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/gamma_pop_densplot_RWcohort", rw_out_tag, ".png"), width = 14, height = 10)

gelman_conv <- cluster_fit$gelman_conv
gelman_multivariate_PSRF <- as.data.frame(gelman_conv$psrf[, 1])
colnames(gelman_multivariate_PSRF) <- "PSRF"
gelman_multivariate_PSRF$param <- c("mu_alpha0", "mu_gamma0",
                                     "overdispersion_parameter", "sd_rw_alpha_pop",
                                     "sd_rw_gamma_pop", "sd_study")
save_psrf_table(cbind.data.frame(parameter = gelman_multivariate_PSRF$param,
                                 PSRF = gelman_multivariate_PSRF$PSRF),
                paste0("./figures_sep22_2026/continuous_cohort_figs/multivariate_PSRF_RW_cohort", rw_out_tag, ".png"))

gelman_par <- cluster_fit$gelman_par
gelman_univariate_PSRF <- as.data.frame(gelman_par$psrf[, 1])
colnames(gelman_univariate_PSRF) <- "PSRF"
gelman_univariate_PSRF$param <- c("log_alpha0_H", "log_alpha0_AIAN", "log_alpha0_API",
                                   "log_alpha0_NHB", "log_alpha0_NHW",
                                   "log_gamma0_H", "log_gamma0_AIAN", "log_gamma0_API",
                                   "log_gamma0_NHB", "log_gamma0_NHW")
save_psrf_table(cbind.data.frame(parameter = gelman_univariate_PSRF$param,
                                 PSRF = gelman_univariate_PSRF$PSRF),
                paste0("./figures_sep22_2026/continuous_cohort_figs/univariate_PSRFF_RW_cohort", rw_out_tag, ".png"))


eff_n_selected <- compute_effective_n(paste0("./results_cluster_sensitivity_analysis/fit_RW_", rw_fit_tag, "continuous_bc.rds"),
                                       data_model_filtered$total)
write_csv(eff_n_selected, paste0("results/effective_sample_size_RWcohort", rw_out_tag, ".csv"))
cat(sprintf(
  "Effective sample size (RW_continuous_bc, overdispersion-adjusted): %.0f (95%% CrI: %.0f, %.0f); raw total_pop = %.0f\n",
  eff_n_selected$eff_n_mean, eff_n_selected$eff_n_LB, eff_n_selected$eff_n_UB, eff_n_selected$total_pop
))

# *****************************************************************************************************
# *****************************************************************************************************
#### MODEL 3: RW cohort, expanded (adds race-specific deviations from the population RW trend) ####
# *****************************************************************************************************
# *****************************************************************************************************

cluster_fit <- readRDS(paste0("./results_cluster_sensitivity_analysis/fit_RW_expanded_", rw_fit_tag, "continuous_bc.rds"))
mcmc_all    <- cluster_fit$mcmc_all
mat_all     <- as.matrix(mcmc_all)

n_iter     <- nrow(mat_all)
cohort_lookup <- make_cohort_lookup()
n_cohorts     <- nrow(cohort_lookup)

# a readable subset of birth cohort years for the by-cohort comparison plots below
# (avoids overplotting when there are many distinct cohort years)
#cohort_years_to_plot <- cohort_lookup$birth_cohort[round(seq(1, n_cohorts, length.out = 6))]

par_alpha0_draws <- t(sapply(1:n_races, function(j) get_vec(mat_all, "par_alpha0", j)))
par_gamma0_draws <- t(sapply(1:n_races, function(j) get_vec(mat_all, "par_gamma0", j)))

rw_alpha_pop_draws <- t(sapply(1:n_cohorts, function(t) get_vec(mat_all, "rw_alpha_pop", t)))
rw_gamma_pop_draws <- t(sapply(1:n_cohorts, function(t) get_vec(mat_all, "rw_gamma_pop", t)))

rw_alpha_dev_draws <- array(NA, dim = c(n_races, n_cohorts, n_iter))
rw_gamma_dev_draws <- array(NA, dim = c(n_races, n_cohorts, n_iter))
for (j in 1:n_races) {
  for (t in 1:n_cohorts) {
    rw_alpha_dev_draws[j, t, ] <- mat_all[, paste0("rw_alpha_dev[", j, ",", t, "]")]
    rw_gamma_dev_draws[j, t, ] <- mat_all[, paste0("rw_gamma_dev[", j, ",", t, "]")]
  }
}

n_studies <- length(unique(data_model_filtered$study_id))
study_re_draws <- t(sapply(1:n_studies, function(s) get_vec(mat_all, "study_re", s)))
theta_draws        <- mat_all[, "theta"]
sd_rw_alpha_draws  <- mat_all[, "sd_rw_alpha_pop"]
sd_rw_gamma_draws  <- mat_all[, "sd_rw_gamma_pop"]

par_fn_rw_exp <- function(r, t) {
  list(
    a_vec = ilogit(par_alpha0_draws[r, ] + rw_alpha_pop_draws[t, ] + rw_alpha_dev_draws[r, t, ]),
    g_vec = exp(   par_gamma0_draws[r, ] + rw_gamma_pop_draws[t, ] + rw_gamma_dev_draws[r, t, ])
  )
}

pred_rw_exp <- run_predictions(par_fn_rw_exp, cohort_lookup, "cohort_idx", race_index, race_labels, ages_vec)
df_p_hat_race_final  <- pred_rw_exp$prev
df_hp_hat_race_final <- pred_rw_exp$foi

for (r in race_index) {
  dt_prev <- df_p_hat_race_final  %>% filter(race == race_labels[r])
  dt_foi  <- df_hp_hat_race_final %>% filter(race == race_labels[r])
  save(dt_prev, dt_foi,
       file = paste0("results/rw_cohort_dt_prev_hat_race_", r, "_",
                      race_short_labels[r], paste0("_RWcohort_expanded", rw_out_tag, ".Rdata")))
}

## cohort trajectory figure (rw_alpha / rw_gamma by year, with race deviations)
rw_full_summary <- expand_grid(race_idx = race_index, cohort_idx = 1:n_cohorts) %>%
  rowwise() %>%
  mutate(
    alpha_mean = mean(rw_alpha_pop_draws[cohort_idx, ] + rw_alpha_dev_draws[race_idx, cohort_idx, ]),
    alpha_lo   = quantile(rw_alpha_pop_draws[cohort_idx, ] + rw_alpha_dev_draws[race_idx, cohort_idx, ], 0.025),
    alpha_hi   = quantile(rw_alpha_pop_draws[cohort_idx, ] + rw_alpha_dev_draws[race_idx, cohort_idx, ], 0.975),
    gamma_mean = mean(rw_gamma_pop_draws[cohort_idx, ] + rw_gamma_dev_draws[race_idx, cohort_idx, ]),
    gamma_lo   = quantile(rw_gamma_pop_draws[cohort_idx, ] + rw_gamma_dev_draws[race_idx, cohort_idx, ], 0.025),
    gamma_hi   = quantile(rw_gamma_pop_draws[cohort_idx, ] + rw_gamma_dev_draws[race_idx, cohort_idx, ], 0.975)
  ) %>%
  ungroup() %>%
  mutate(race = race_labels[race_idx]) %>%
  left_join(cohort_lookup, by = "cohort_idx")

p_rw_alpha <- ggplot(rw_full_summary, aes(x = birth_cohort, y = alpha_mean, color = race, fill = race)) +
  geom_ribbon(aes(ymin = alpha_lo, ymax = alpha_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  # ticks mark birth years with data; years between them are RW-interpolated
  geom_rug(data = filter(cohort_lookup, observed), aes(x = birth_cohort),
           inherit.aes = FALSE, sides = "b", alpha = 0.6) +
  labs(x = "Birth cohort year", y = "Cohort effect on alpha (logit scale)",
       title = "Estimated cohort trajectory: max-prevalence (alpha) component, by race") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

p_rw_gamma <- ggplot(rw_full_summary, aes(x = birth_cohort, y = gamma_mean, color = race, fill = race)) +
  geom_ribbon(aes(ymin = gamma_lo, ymax = gamma_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  # ticks mark birth years with data; years between them are RW-interpolated
  geom_rug(data = filter(cohort_lookup, observed), aes(x = birth_cohort),
           inherit.aes = FALSE, sides = "b", alpha = 0.6) +
  labs(x = "Birth cohort year", y = "Cohort effect on gamma (log scale)",
       title = "Estimated cohort trajectory: rate (gamma) component, by race") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/rw_alpha_cohort_trajectory", rw_out_tag, ".png"), p_rw_alpha, width = 10, height = 6)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/rw_gamma_cohort_trajectory", rw_out_tag, ".png"), p_rw_gamma, width = 10, height = 6)

df_foi_curves_subset <- df_hp_hat_race_final %>%
  filter( age > 0, Mean < 1)

ggplot(df_foi_curves_subset,
       aes(x = age, y = Mean, color = birth_cohort, group = birth_cohort)) +
  geom_line(linewidth = 1) +
  facet_wrap(~race, ncol = 3) +
  scale_color_viridis_c(name = "Birth cohort\nyear") +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  xlab("Age") + ylab("Predicted force of infection") +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "bold",size=20))

ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/predicted_foi_by_cohort_year_RWcohort_expanded", rw_out_tag, ".png"), width = 12, height = 8)

df_p_curves_subset <- df_p_hat_race_final %>% filter(age > 0)
ggplot(df_p_curves_subset,
       aes(x = age, y = Mean, color = birth_cohort, group = birth_cohort)) +
  geom_line(linewidth = 1) +
  facet_wrap(~race, ncol = 3) +
  scale_color_viridis_c(name = "Birth cohort\nyear") +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  xlab("Age") + ylab("Predicted seroprevalence") +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "bold",size=20))

ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/predicted_prev_by_cohort_year_RWcohort_expanded", rw_out_tag, ".png"), width = 12, height = 8)


max_prev_by_race_cohort <- df_p_hat_race_final %>%
  filter(age == fin_age)#, birth_cohort %in% cohort_years_to_plot)

ggplot(max_prev_by_race_cohort,
       aes(x = race, y = Mean, color = factor(birth_cohort),
           ymin = LB, ymax = UB, group = factor(birth_cohort))) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.7) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
  labs(x = "", y = paste0("Max seroprevalence (age ", fin_age, ")"),
       color = "Birth cohort\nyear") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/max_prev_by_race_cohort_RWcohort_expanded", rw_out_tag, ".png"), width = 11, height = 8)

plot_max_prev_rw_exp <- function(df, ref_age) {
  df <- df%>%
    mutate(Mean=Mean*100)%>%
    mutate(LB=LB*100)%>%
    mutate(UB=UB*100)
  ggplot(df,
         aes(x = race, y = Mean, color = factor(birth_cohort),
             ymin = LB, ymax = UB, group = factor(birth_cohort))) +
    geom_pointrange(position = position_dodge(width = 0.6), size = 0.7) +
    scale_y_continuous(breaks = seq(0,100, by = 10), limits = c(0, 100)) +
    labs(x = "", y = paste0("Seroprevalence (age ", ref_age, ")"),
         color = "Birth cohort\nyear") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}

max_prev_by_race_cohort_list <- lapply(ref_ages, function(a)
  make_max_prev_by_race_cohort(df_p_hat_race_final, a))
names(max_prev_by_race_cohort_list) <- paste0("age_", ref_ages)
max_prev_by_race_cohort <- max_prev_by_race_cohort_list[[paste0("age_", fin_age)]]

for (ref_age in ref_ages) {
  df_this <- max_prev_by_race_cohort_list[[paste0("age_", ref_age)]]
  p_this  <- plot_max_prev_rw_exp(df_this, ref_age)

  out_file <- if (ref_age == fin_age) {
    paste0("./figures_sep22_2026/continuous_cohort_figs/max_prev_by_race_cohort_RWcohort_expanded", rw_out_tag, ".png")
  } else {
    paste0("./figures_sep22_2026/continuous_cohort_figs//max_prev_by_race_cohort_age", ref_age, "_RWcohort_expanded", rw_out_tag, ".png")
  }
  ggsave(out_file, p_this, width = 11, height = 8)
}

max_prev_by_race_cohort_all_ages <- bind_rows(max_prev_by_race_cohort_list)
save(max_prev_by_race_cohort_all_ages,
     file = paste0("results/rw_cohort_expanded_max_prev_age5_10_90", rw_out_tag, ".Rdata"))

## median age of infection
gamma_fn_rw_exp <- function(r, t, iter) {
  exp(par_gamma0_draws[r, iter] + rw_gamma_pop_draws[t, iter] + rw_gamma_dev_draws[r, t, iter])
}
df_aoi_final <- run_aoi(gamma_fn_rw_exp, cohort_lookup, "cohort_idx", race_index, race_labels, n_iter)

for (r in race_index) {
  dt_aoi <- df_aoi_final %>% filter(race == race_labels[r])
  save(dt_aoi, file = paste0("results/rw_cohort_dt_aoi_race_", r, "_",
                              race_short_labels[r], paste0("_RWcohort_expanded", rw_out_tag, ".Rdata")))
}

save(df_aoi_final, file = paste0("results/rw_cohort_expanded_df_aoi_all", rw_out_tag, ".Rdata"))

aoi_by_race_cohort <- df_aoi_final %>%
  filter(Mean < 100)
  #filter(birth_cohort %in% cohort_years_to_plot)

ggplot(aoi_by_race_cohort,
       aes(x = race, y = Mean, color = factor(birth_cohort),
           ymin = LB, ymax = UB, group = factor(birth_cohort))) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.7) +
  labs(x = "", y = "Median age of infection (years)", color = "Birth cohort\nyear") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/median_age_infection_by_race_cohort_RWcohort_expanded", rw_out_tag, ".png"), width = 11, height = 8)

## key parameter summary
param_summary <- tibble(
  parameter = c("sd_rw_alpha_pop", "sd_rw_gamma_pop", "theta"),
  Mean = c(mean(sd_rw_alpha_draws), mean(sd_rw_gamma_draws), mean(theta_draws)),
  LB   = c(quantile(sd_rw_alpha_draws, 0.025), quantile(sd_rw_gamma_draws, 0.025), quantile(theta_draws, 0.025)),
  UB   = c(quantile(sd_rw_alpha_draws, 0.975), quantile(sd_rw_gamma_draws, 0.975), quantile(theta_draws, 0.975))
)
print(param_summary)
write_csv(param_summary, paste0("results/rw_cohort_key_param_summary", rw_out_tag, ".csv"))

## convergence plots/tables
m_alpha0 <- mat_all[, paste0("par_alpha0[", 1:n_races, "]"), drop = FALSE]
m_gamma0 <- mat_all[, paste0("par_gamma0[", 1:n_races, "]"), drop = FALSE]
alpha0_labels <- paste0("alpha intercept ", race_labels)
gamma0_labels <- paste0("gamma intercept ", race_labels)

generate_traceplot_function(m_alpha0, alpha0_labels, ncol_facet = 1)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/alpha0_traceplot_RWcohort_expanded", rw_out_tag, ".png"), width = 12, height = 8)
generate_traceplot_function(m_gamma0, gamma0_labels, ncol_facet = 1)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/gamma0_traceplot_RWcohort_expanded", rw_out_tag, ".png"), width = 12, height = 8)
generate_densplot_function(m_alpha0, alpha0_labels, ncol_facet = 1)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/alpha0_densplot_RWcohort_expanded", rw_out_tag, ".png"), width = 12, height = 8)
generate_densplot_function(m_gamma0, gamma0_labels, ncol_facet = 1)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/gamma0_densplot_RWcohort_expanded", rw_out_tag, ".png"), width = 12, height = 8)

ref_cohort_bins  <- round(seq(1, n_cohorts, length.out = 6))
ref_cohort_years <- cohort_lookup$birth_cohort[ref_cohort_bins]

m_alpha_pop <- mat_all[, paste0("rw_alpha_pop[", ref_cohort_bins, "]"), drop = FALSE]
m_gamma_pop <- mat_all[, paste0("rw_gamma_pop[", ref_cohort_bins, "]"), drop = FALSE]
alpha_pop_labels <- paste0("alpha_pop, cohort=", ref_cohort_years)
gamma_pop_labels <- paste0("gamma_pop, cohort=", ref_cohort_years)

generate_traceplot_function(m_alpha_pop, alpha_pop_labels, ncol_facet = 2)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/alpha_pop_traceplot_RWcohort_expanded", rw_out_tag, ".png"), width = 14, height = 10)
generate_traceplot_function(m_gamma_pop, gamma_pop_labels, ncol_facet = 2)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/gamma_pop_traceplot_RWcohort_expanded", rw_out_tag, ".png"), width = 14, height = 10)
generate_densplot_function(m_alpha_pop, alpha_pop_labels, ncol_facet = 2)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/alpha_pop_densplot_RWcohort_expanded", rw_out_tag, ".png"), width = 14, height = 10)
generate_densplot_function(m_gamma_pop, gamma_pop_labels, ncol_facet = 2)
ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/gamma_pop_densplot_RWcohort_expanded", rw_out_tag, ".png"), width = 14, height = 10)

for (r in 1:n_races) {
  m_alpha_dev_r <- mat_all[, paste0("rw_alpha_dev[", r, ",", ref_cohort_bins, "]"), drop = FALSE]
  m_gamma_dev_r <- mat_all[, paste0("rw_gamma_dev[", r, ",", ref_cohort_bins, "]"), drop = FALSE]

  alpha_dev_labels <- paste0(race_labels[r], " alpha_dev, cohort=", ref_cohort_years)
  gamma_dev_labels <- paste0(race_labels[r], " gamma_dev, cohort=", ref_cohort_years)

  race_short_r <- race_short_labels[r]  # e.g. "AIAN", not "AI/AN" -- a literal "/" would break the output path

  generate_traceplot_function(m_alpha_dev_r, alpha_dev_labels, ncol_facet = 2)
  ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/alpha_dev_traceplot_", race_short_r, "_RWcohort_expanded", rw_out_tag, ".png"), width = 14, height = 10)
  generate_traceplot_function(m_gamma_dev_r, gamma_dev_labels, ncol_facet = 2)
  ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/gamma_dev_traceplot_", race_short_r, "_RWcohort_expanded", rw_out_tag, ".png"), width = 14, height = 10)
  generate_densplot_function(m_alpha_dev_r, alpha_dev_labels, ncol_facet = 2)
  ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/alpha_dev_densplot_", race_short_r, "_RWcohort_expanded", rw_out_tag, ".png"), width = 14, height = 10)
  generate_densplot_function(m_gamma_dev_r, gamma_dev_labels, ncol_facet = 2)
  ggsave(paste0("./figures_sep22_2026/continuous_cohort_figs/gamma_dev_densplot_", race_short_r, "_RWcohort_expanded", rw_out_tag, ".png"), width = 14, height = 10)
}

gelman_conv <- cluster_fit$gelman_conv
gelman_multivariate_PSRF <- as.data.frame(gelman_conv$psrf[, 1])
colnames(gelman_multivariate_PSRF) <- "PSRF"
gelman_multivariate_PSRF$param <- c("mu_alpha0", "mu_gamma0",
                                     "overdispersion_parameter", "sd_rw_alpha_pop",
                                     "sd_rw_gamma_pop", "sd_rw_dev_alpha",
                                     "sd_rw_dev_gamma", "sd_study")
save_psrf_table(cbind.data.frame(parameter = gelman_multivariate_PSRF$param,
                                 PSRF = gelman_multivariate_PSRF$PSRF),
                paste0("./figures_sep22_2026/continuous_cohort_figs/multivariate_PSRF_RW_cohort_expanded", rw_out_tag, ".png"))

gelman_par <- cluster_fit$gelman_par
gelman_univariate_PSRF <- as.data.frame(gelman_par$psrf[, 1])
colnames(gelman_univariate_PSRF) <- "PSRF"
gelman_univariate_PSRF$param <- c("log_alpha0_H", "log_alpha0_AIAN", "log_alpha0_API",
                                   "log_alpha0_NHB", "log_alpha0_NHW",
                                   "log_gamma0_H", "log_gamma0_AIAN", "log_gamma0_API",
                                   "log_gamma0_NHB", "log_gamma0_NHW")
save_psrf_table(cbind.data.frame(parameter = gelman_univariate_PSRF$param,
                                 PSRF = gelman_univariate_PSRF$PSRF),
                paste0("./figures_sep22_2026/continuous_cohort_figs/univariate_PSRFF_RW_cohort_expanded", rw_out_tag, ".png"))


