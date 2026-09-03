# *****************************************************************************
#### Load necessary packages ####
# *****************************************************************************

library(tidyverse)
library(data.table)
library(coda)
library(matrixStats)
library(gt)

# *****************************************************************************
#### Data prep ####
# *****************************************************************************

load("data/df_hp_usa_w_final.Rdata")

data_model <- df_hp_usa_w_final %>%
  arrange(race_f, year) %>%
  filter(race_f != 'Other') %>%
  mutate(study_id_fac = factor(study_id),
         race_fac = factor(race_f),
         time = year - 1988,
         bc_scaled = birth_cohort - min(birth_cohort)) %>%
  mutate(birth_cohort_cat = ifelse(birth_cohort < 1960, 0, 1))

data_model_filtered <- data_model %>%
  filter(study_id != "SawJRhe2004") %>%
  filter(!(study_id == 'GooJImm2003' & subgroup_n == 384.00)) %>%
  filter(!(study_id == 'FelBirt2007' & subgroup_n == 147.00)) %>%
  filter(!(study_id == 'ParClin2000' & subgroup_n == 2080.00)) %>%
  filter(!(study_id == 'ParClin2000' & subgroup_n == 1040.00)) %>%
  filter(!(study_id == 'CorCanc1990' & age == 0.0410)) %>%
  filter(!(study_id == 'CorCanc1990' & age == 0.2500)) %>%
  filter(!(study_id == 'CorCanc1990' & age == 0.7085))

cluster_fit <- readRDS("./cluster_output/fit_base_continuous_bc.rds")
mcmc_all    <- cluster_fit$mcmc_all
mat_all     <- as.matrix(mcmc_all)

ilogit <- function(logodds) 1 / (1 + exp(-logodds))

n_races <- 5
fin_age <- 90
n_iter  <- nrow(mat_all)

race_index  <- 1:5
race_labels <- c("Hispanic", "NH American Indian", "NH Asian", "NH Black", "NH White")
race_short  <- c("H", "AIAN", "NHA", "NHB", "NHW")

dir.create("figures_aug20_2026", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# *****************************************************************************
#### Helper functions ####
# *****************************************************************************

## representative birth-cohort years to evaluate the linear predictor at
## (replaces cohort_idx/n_cohorts entirely -- this model has no cohort bins,
## bc_scaled is a continuous covariate)
bc_range <- range(data_model_filtered$bc_scaled, na.rm = TRUE)
n_bc_points <- 10   # number of representative points across the cohort range

bc_scaled_grid <- seq(bc_range[1], bc_range[2], length.out = n_bc_points)
birth_cohort_grid <- bc_scaled_grid + min(data_model_filtered$birth_cohort)

bc_lookup <- data.frame(
  bc_point_idx  = seq_len(n_bc_points),
  bc_scaled     = bc_scaled_grid,
  birth_cohort  = round(birth_cohort_grid)
)

## process results output
get_vec <- function(mat_all, base_name, idx) {
  mat_all[, paste0(base_name, "[", idx, "]")]
}

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

# -----------------------------------------------------------------------------
# NOTE on study_re in the prediction curves below: predictions are POPULATION
# AVERAGE curves (study_re marginalized to 0), consistent with the earlier
# convention in this project -- study_re represents a study-specific shift,
# not part of the race/cohort-level trend being visualized here.
# -----------------------------------------------------------------------------

# *****************************************************************************
#### Prediction loop -- all races, evaluated at representative bc_scaled points ####
# *****************************************************************************

ages_vec <- 0:fin_age
n_age    <- fin_age + 1

prev_groups <- vector("list", length(race_index))
foi_groups  <- vector("list", length(race_index))

for (r in race_index) {

  prev_mean <- prev_sd <- prev_lb <- prev_ub <- matrix(NA_real_, n_age, n_bc_points)
  foi_mean  <- foi_sd  <- foi_lb  <- foi_ub  <- matrix(NA_real_, n_age, n_bc_points)

  for (t in 1:n_bc_points) {

    bc_val <- bc_scaled_grid[t]

    a_vec <- ilogit(m_alpha0[r, ] + m_alpha1[r, ] * bc_val)
    g_vec <- exp(   m_gamma0[r, ] + m_gamma1[r, ] * bc_val)

    g_age         <- outer(ages_vec, g_vec, FUN = "*")
    exp_neg       <- exp(-g_age)
    one_minus_exp <- 1 - exp_neg

    a_rep <- rep(a_vec, each = n_age)
    g_rep <- rep(g_vec, each = n_age)

    p_raw_mat <- one_minus_exp * a_rep
    p_mat     <- 0.001 + 0.998 * p_raw_mat

    foi_num <- a_rep * g_rep * exp_neg
    foi_den <- 1 - p_mat * one_minus_exp + 1e-10
    foi_mat <- foi_num / foi_den
    foi_mat[1, ] <- NA

    prev_mean[, t] <- rowMeans2(p_raw_mat, na.rm = TRUE)
    prev_sd[, t]   <- rowSds(p_raw_mat,   na.rm = TRUE)
    pq <- rowQuantiles(p_raw_mat, probs = c(0.025, 0.975), na.rm = TRUE)
    prev_lb[, t] <- pq[, 1]; prev_ub[, t] <- pq[, 2]

    foi_mean[, t] <- rowMeans2(foi_mat, na.rm = TRUE)
    foi_sd[, t]   <- rowSds(foi_mat,   na.rm = TRUE)
    fq <- rowQuantiles(foi_mat, probs = c(0.025, 0.975), na.rm = TRUE)
    foi_lb[, t] <- fq[, 1]; foi_ub[, t] <- fq[, 2]
  }

  to_long <- function(mat, stat_name) {
    as.data.table(reshape2::melt(mat, varnames = c("age_idx", "bc_point_idx"),
                                  value.name = stat_name))
  }

  dt_prev <- to_long(prev_mean, "Mean")[to_long(prev_sd, "SD"), on = c("age_idx", "bc_point_idx")
           ][to_long(prev_lb, "LB"), on = c("age_idx", "bc_point_idx")
           ][to_long(prev_ub, "UB"), on = c("age_idx", "bc_point_idx")]
  dt_prev[, `:=`(age = age_idx - 1, race = race_labels[r], type = "Prevalence")]
  dt_prev <- merge(dt_prev, bc_lookup, by = "bc_point_idx", all.x = TRUE)

  dt_foi <- to_long(foi_mean, "Mean")[to_long(foi_sd, "SD"), on = c("age_idx", "bc_point_idx")
          ][to_long(foi_lb, "LB"), on = c("age_idx", "bc_point_idx")
          ][to_long(foi_ub, "UB"), on = c("age_idx", "bc_point_idx")]
  dt_foi[, `:=`(age = age_idx - 1, race = race_labels[r], type = "Force of Infection")]
  dt_foi <- merge(dt_foi, bc_lookup, by = "bc_point_idx", all.x = TRUE)

  prev_groups[[r]] <- dt_prev
  foi_groups[[r]]  <- dt_foi

  save(dt_prev, dt_foi,
       file = paste0("results/base_continuous_dt_prev_hat_race_", r, "_",
                      race_short[r], "_BaseContinuous.Rdata"))

  cat("Finished predictions for race:", race_labels[r], "\n")
}

df_p_hat_race_final  <- bind_rows(prev_groups)
df_hp_hat_race_final <- bind_rows(foi_groups)

# *****************************************************************************************************
#### FOI and prevalence curves across ages by race and birth cohort year ####
# *****************************************************************************************************

cohort_years_to_plot <- bc_lookup$birth_cohort  # all n_bc_points representative years

df_foi_curves_subset <- df_hp_hat_race_final %>%
  filter(age > 0, Mean < 1)

ggplot(df_foi_curves_subset,
       aes(x = age, y = Mean, color = birth_cohort, group = birth_cohort)) +
  geom_line(linewidth = 1) +
  facet_wrap(~race, ncol = 3) +
  scale_color_viridis_c(name = "Birth cohort\nyear") +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  xlab("Age") + ylab("Predicted force of infection") +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "bold"))

ggsave("./figures_aug20_2026/predicted_foi_by_cohort_year_BaseContinuous.png",
       width = 12, height = 8)

df_p_curves_subset <- df_p_hat_race_final %>%
  filter(age > 0)

ggplot(df_p_curves_subset,
       aes(x = age, y = Mean, color = birth_cohort, group = birth_cohort)) +
  geom_line(linewidth = 1) +
  facet_wrap(~race, ncol = 3) +
  scale_color_viridis_c(name = "Birth cohort\nyear") +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  xlab("Age") + ylab("Predicted seroprevalence") +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "bold"))

ggsave("./figures_aug20_2026/predicted_prev_by_cohort_year_BaseContinuous.png",
       width = 12, height = 8)

# *****************************************************************************
#### Max prevalence at multiple ages ####
# *****************************************************************************
make_max_prev_by_race_cohort <- function(ref_age) {
  df_p_hat_race_final %>%
    filter(age == ref_age) %>%
    mutate(ref_age = ref_age)
}

plot_max_prev_by_race_cohort <- function(df, ref_age) {
  ggplot(df,
         aes(x = race, y = Mean, color = factor(birth_cohort),
             ymin = LB, ymax = UB, group = factor(birth_cohort))) +
    geom_pointrange(position = position_dodge(width = 0.6), size = 0.7) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
    labs(x = "", y = paste0("Seroprevalence (age ", ref_age, ")"),
         color = "Birth cohort\nyear") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
}

ref_ages <- c(5, 10, fin_age)

max_prev_by_race_cohort_list <- lapply(ref_ages, make_max_prev_by_race_cohort)
names(max_prev_by_race_cohort_list) <- paste0("age_", ref_ages)
max_prev_by_race_cohort <- max_prev_by_race_cohort_list[[paste0("age_", fin_age)]]

for (ref_age in ref_ages) {
  df_this <- max_prev_by_race_cohort_list[[paste0("age_", ref_age)]]
  p_this  <- plot_max_prev_by_race_cohort(df_this, ref_age)

  out_file <- if (ref_age == fin_age) {
    "./figures_aug20_2026/max_prev_by_race_cohort_BaseContinuous.png"
  } else {
    paste0("./figures_aug20_2026/max_prev_by_race_cohort_age", ref_age, "_BaseContinuous.png")
  }

  ggsave(out_file, p_this, width = 11, height = 8)
}

max_prev_by_race_cohort_all_ages <- bind_rows(max_prev_by_race_cohort_list)
save(max_prev_by_race_cohort_all_ages,
     file = "results/base_continuous_max_prev_age5_10_90.Rdata")

# *****************************************************************************
#### Computing median age of infection ####
# *****************************************************************************

aoi_groups <- list()

for (r in race_index) {

  aoi_mat <- matrix(NA, nrow = n_iter, ncol = n_bc_points)

  for (t in 1:n_bc_points) {
    bc_val <- bc_scaled_grid[t]
    for (iter in 1:n_iter) {
      g <- exp(m_gamma0[r, iter] + m_gamma1[r, iter] * bc_val)
      aoi_mat[iter, t] <- log(2) / g
    }
  }

  dt_aoi <- data.frame(
    bc_point_idx = 1:n_bc_points,
    Mean = apply(aoi_mat, 2, mean),
    LB   = apply(aoi_mat, 2, quantile, 0.025),
    UB   = apply(aoi_mat, 2, quantile, 0.975)
  ) %>%
    left_join(bc_lookup, by = "bc_point_idx") %>%
    mutate(race = race_labels[r], type = "Median age of Infection")

  aoi_groups[[r]] <- dt_aoi

  save(dt_aoi,
       file = paste0("results/base_continuous_dt_aoi_race_", r, "_",
                      race_short[r], "_BaseContinuous.Rdata"))

  rm(aoi_mat, dt_aoi)
  cat("Finished median age of infection for race:", race_labels[r], "\n")
}

df_aoi_final <- bind_rows(aoi_groups)

write_csv(df_aoi_final, "results/median_age_inf_BaseContinuous.csv")

# *****************************************************************************
#### Visualizations of model convergence ####
# *****************************************************************************

alpha0_labels <- paste0("alpha intercept ", race_labels)
alpha1_labels <- paste0("alpha slope ", race_labels)
gamma0_labels <- paste0("gamma intercept ", race_labels)
gamma1_labels <- paste0("gamma slope ", race_labels)

## generate traceplots
generate_traceplot_function <- function(draws, labels, ncol_facet = 1) {
  draws_df <- as.data.frame(draws)
  colnames(draws_df) <- labels
  draws_df$iteration <- seq_len(nrow(draws_df))

  draws_df_long <- pivot_longer(draws_df,
                                 cols = -iteration,
                                 names_to = "parameter",
                                 values_to = "value")

  ggplot(draws_df_long, aes(x = iteration, y = value, color = parameter)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ parameter, scales = "free_y", ncol = ncol_facet) +
    theme_minimal() +
    labs(x = "Iteration", y = "Estimate") +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
          strip.text = element_text(size = 10))
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
          axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
          strip.text = element_text(size = 10))
}

# note: m_alpha0 etc. are [race x iter] here (transposed above) -- transpose
# back to [iter x race] for these plotting functions, which expect draws in rows
generate_traceplot_function(t(m_alpha0), alpha0_labels, ncol_facet = 1)
ggsave("./figures_aug20_2026/alpha0_traceplot_BaseContinuous.png", width = 12, height = 8)
generate_traceplot_function(t(m_alpha1), alpha1_labels, ncol_facet = 1)
ggsave("./figures_aug20_2026/alpha1_traceplot_BaseContinuous.png", width = 12, height = 8)
generate_traceplot_function(t(m_gamma0), gamma0_labels, ncol_facet = 1)
ggsave("./figures_aug20_2026/gamma0_traceplot_BaseContinuous.png", width = 12, height = 8)
generate_traceplot_function(t(m_gamma1), gamma1_labels, ncol_facet = 1)
ggsave("./figures_aug20_2026/gamma1_traceplot_BaseContinuous.png", width = 12, height = 8)

generate_densplot_function(t(m_alpha0), alpha0_labels, ncol_facet = 1)
ggsave("./figures_aug20_2026/alpha0_densplot_BaseContinuous.png", width = 12, height = 8)
generate_densplot_function(t(m_alpha1), alpha1_labels, ncol_facet = 1)
ggsave("./figures_aug20_2026/alpha1_densplot_BaseContinuous.png", width = 12, height = 8)
generate_densplot_function(t(m_gamma0), gamma0_labels, ncol_facet = 1)
ggsave("./figures_aug20_2026/gamma0_densplot_BaseContinuous.png", width = 12, height = 8)
generate_densplot_function(t(m_gamma1), gamma1_labels, ncol_facet = 1)
ggsave("./figures_aug20_2026/gamma1_densplot_BaseContinuous.png", width = 12, height = 8)

## generate tables for univariate and multivariate convergence statistics
## (this model's ForMCMCConv block: mub[1:4], var_b[1:4], cov_b[1:6],
##  theta, sd_study = 16 parameters -- matches model_base_continuous_bc)

gelman_conv <- cluster_fit$gelman_conv

gelman_multivariate_PSRF <- as.data.frame(gelman_conv$psrf[, 1])
colnames(gelman_multivariate_PSRF) <- "PSRF"
gelman_multivariate_PSRF$param <- c(
  "mu_intercept_alpha", "mu_slope_alpha", "mu_intercept_gamma", "mu_slope_gamma",
  "var_intercept_alpha", "var_slope_alpha", "var_intercept_gamma", "var_slope_gamma",
  "cov_slope_intercept_alpha", "cov_intercept_gamma_intercept_alpha",
  "cov_intercept_gamma_slope_alpha", "cov_slope_gamma_intercept_alpha",
  "cov_slope_gamma_slope_alpha", "cov_slope_intercept_gamma",
  "overdispersion_parameter", "sd_study"
)

gelman_multivariate_PSRF_table <- gt(cbind.data.frame(
  parameter = gelman_multivariate_PSRF$param,
  PSRF      = gelman_multivariate_PSRF$PSRF
))
gelman_multivariate_PSRF_table
gtsave(gelman_multivariate_PSRF_table, "./figures_aug20_2026/multivariate_PSRF_BaseContinuous.png")

gelman_par <- cluster_fit$gelman_par
gelman_univariate_PSRF <- as.data.frame(gelman_par$psrf[, 1])
colnames(gelman_univariate_PSRF) <- "PSRF"
gelman_univariate_PSRF$param <- c(
  "log_odds_alpha_H", "log_odds_alpha_AIAN", "log_odds_alpha_NHA",
  "log_odds_alpha_NHB", "log_odds_alpha_NHW",
  "delta_log_odds_alpha_H", "delta_log_odds_alpha_AIAN", "delta_log_odds_alpha_NHA",
  "delta_log_odds_alpha_NHB", "delta_log_odds_alpha_NHW",
  "log_gamma_H", "log_gamma_AIAN", "log_gamma_NHA",
  "log_gamma_NHB", "log_gamma_NHW",
  "delta_log_gamma_H", "delta_log_gamma_AIAN", "delta_log_gamma_NHA",
  "delta_log_gamma_NHB", "delta_log_gamma_NHW"
)

gelman_univariate_PSRF_table <- gt(cbind.data.frame(
  parameter = gelman_univariate_PSRF$param,
  PSRF      = gelman_univariate_PSRF$PSRF
))
gelman_univariate_PSRF_table
gtsave(gelman_univariate_PSRF_table, "./figures_aug20_2026/univariate_PSRF_BaseContinuous.png")

cat("\nDone. All figures/tables saved to figures_aug20_2026/.\n")