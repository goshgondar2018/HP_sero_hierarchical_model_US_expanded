library(tidyverse)
library(data.table)
library(coda)
library(matrixStats)

# load data
load("data/df_hp_usa_w_final.Rdata")

# sort data by year of the study and study id
data_model <- df_hp_usa_w_final %>%
  arrange(race_f, year)

data_model <- data_model %>%
  filter(race_f != 'Other')
# %>% filter(age > 0.5)

data_model <- data_model %>%
  mutate(study_id_fac = factor(study_id),
         race_fac = factor(race_f),
         time = year - 1988,
         bc_scaled = birth_cohort - min(birth_cohort))%>%
  mutate(birth_cohort_cat=ifelse(birth_cohort<1960,0,1))

# remove aberrant measures for NHB population
#data_model_filtered=data_model[c(1:171,174:nrow(data_model)),]

# remove study duplicates
data_model_filtered=data_model%>%
  filter(study_id!="SawJRhe2004")%>%
  filter(!(study_id=='GooJImm2003'&subgroup_n==384.00))%>%
  filter(!(study_id=='FelBirt2007'&subgroup_n==147.00))%>%
  filter(!(study_id=='ParClin2000'&subgroup_n==2080.00))%>%
  filter(!(study_id=='ParClin2000'&subgroup_n==1040.00))%>%
  filter(!(study_id=='CorCanc1990'&age==0.0410))%>%
  filter(!(study_id=='CorCanc1990'&age==0.2500))%>%
  filter(!(study_id=='CorCanc1990'&age==0.7085))

cluster_fit <- readRDS("./cluster_output/fit_RW_expanded_continuous_bc.rds")
mcmc_all    <- cluster_fit$mcmc_all
mat_all     <- as.matrix(mcmc_all)

ilogit <- function(logodds) 1 / (1 + exp(-logodds))

fin_age <- 90

race_index  <- 1:5
race_labels <- c("Hispanic", "NH American Indian", "NH Asian", "NH Black", "NH White")
race_short  <- c("H", "AIAN", "NHA", "NHB", "NHW")
cohort_idx <- as.numeric(as.factor(round(data_model_filtered$birth_cohort)))
N <- nrow(data_model_filtered)
inds <- as.numeric(data_model_filtered$race_fac)
n_races <- length(unique(data_model_filtered$race_fac))
hp <- round(data_model_filtered$HP)
age <- data_model_filtered$age  
bc_scaled <- data_model_filtered$bc_scaled
birth_cohort_cat=data_model_filtered$birth_cohort_cat
n <- round(data_model_filtered$subgroup_n)
n_cohorts <- nrow(cohort_lookup)
study_id_raw <- data_model_filtered$study_id
study_id     <- as.numeric(as.factor(study_id_raw))

# -----------------------------------------------------------------------------
# 1. Cohort bin <-> actual birth-cohort-year lookup
#    (rw_alpha/rw_gamma are indexed 1:n_cohorts, not by calendar year)
# -----------------------------------------------------------------------------
cohort_lookup <- data.frame(
  cohort_idx   = cohort_idx,
  birth_cohort = round(data_model_filtered$birth_cohort)
) %>%
  distinct() %>%
  arrange(cohort_idx)


# -----------------------------------------------------------------------------
# 2. Extract posterior draws from fit_model (jags.samples()-style mcarray)
#
#    Expected dims:
#      fit_model[["par_alpha0"]] : [n_races, iter_per_chain, n_chains]
#      fit_model[["par_gamma0"]] : [n_races, iter_per_chain, n_chains]
#      fit_model[["rw_alpha"]]   : [n_cohorts, iter_per_chain, n_chains]
#      fit_model[["rw_gamma"]]   : [n_cohorts, iter_per_chain, n_chains]
#      fit_model[["theta"]]      : [1, iter_per_chain, n_chains]
#
#    IMPORTANT: confirm with dim(fit_model[["par_alpha0"]]) etc. before
#    trusting this blindly -- if your actual object came from runjags /
#    as.mcmc.list() instead of jags.samples(), the structure will be a flat
#    matrix of named columns instead, and this extraction section needs to
#    be rewritten (ask if that's the case).
# -----------------------------------------------------------------------------
n_iter <- nrow(mat_all)

get_vec <- function(mat_all, base_name, idx) {
  mat_all[, paste0(base_name, "[", idx, "]")]
}

# par_alpha0[1..n_races], par_gamma0[1..n_races]: [n_iter x n_races]
par_alpha0_draws <- sapply(1:n_races, function(j) get_vec(mat_all, "par_alpha0", j))
par_gamma0_draws <- sapply(1:n_races, function(j) get_vec(mat_all, "par_gamma0", j))
# transpose to match your original [race, iter] orientation used below
par_alpha0_draws <- t(par_alpha0_draws)
par_gamma0_draws <- t(par_gamma0_draws)

# rw_alpha_pop[1..n_cohorts], rw_gamma_pop[1..n_cohorts]: [n_cohorts x n_iter]
rw_alpha_pop_draws <- t(sapply(1:n_cohorts, function(t) get_vec(mat_all, "rw_alpha_pop", t)))
rw_gamma_pop_draws <- t(sapply(1:n_cohorts, function(t) get_vec(mat_all, "rw_gamma_pop", t)))

# rw_alpha_dev[race,cohort], rw_gamma_dev[race,cohort]: [n_races, n_cohorts, n_iter]
rw_alpha_dev_draws <- array(NA, dim = c(n_races, n_cohorts, n_iter))
rw_gamma_dev_draws <- array(NA, dim = c(n_races, n_cohorts, n_iter))
for (j in 1:n_races) {
  for (t in 1:n_cohorts) {
    rw_alpha_dev_draws[j, t, ] <- mat_all[, paste0("rw_alpha_dev[", j, ",", t, "]")]
    rw_gamma_dev_draws[j, t, ] <- mat_all[, paste0("rw_gamma_dev[", j, ",", t, "]")]
  }
}

# study_re[1..n_studies]: [n_studies x n_iter]
n_studies <- length(unique(study_id))
study_re_draws <- t(sapply(1:n_studies, function(s) get_vec(mat_all, "study_re", s)))

theta_draws <- mat_all[, "theta"]
sd_rw_alpha_draws <- mat_all[, "sd_rw_alpha_pop"]   # note: name matches model's actual node
sd_rw_gamma_draws <- mat_all[, "sd_rw_gamma_pop"]

# -----------------------------------------------------------------------------
# 3. Prediction loop: race x cohort_bin x age
# -----------------------------------------------------------------------------

ages_vec <- 0:fin_age
n_age    <- fin_age + 1

prev_groups <- vector("list", length(race_index))
foi_groups  <- vector("list", length(race_index))

for (r in race_index) {

  prev_mean <- prev_sd <- prev_lb <- prev_ub <- matrix(NA_real_, n_age, n_cohorts)
  foi_mean  <- foi_sd  <- foi_lb  <- foi_ub  <- matrix(NA_real_, n_age, n_cohorts)

  for (t in 1:n_cohorts) {

    a_vec <- ilogit(par_alpha0_draws[r, ] + rw_alpha_pop_draws[t, ] +
                     rw_alpha_dev_draws[r, t, ])
    g_vec <- exp(   par_gamma0_draws[r, ] + rw_gamma_pop_draws[t, ] +
                     rw_gamma_dev_draws[r, t, ])

    g_age         <- outer(ages_vec, g_vec, FUN = "*")
    exp_neg       <- exp(-g_age)
    one_minus_exp <- 1 - exp_neg

    a_rep <- rep(a_vec, each = n_age)   # broadcasts a_vec down each column
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

  # everything below is unchanged — prev_mean/prev_sd/... and
  # foi_mean/foi_sd/... are already exactly the matrices your
  # to_long()/merge()/save() block expects
  to_long <- function(mat, stat_name) {
    as.data.table(reshape2::melt(mat, varnames = c("age_idx", "cohort_idx"),
                                  value.name = stat_name))
  }

  dt_prev <- to_long(prev_mean, "Mean")[to_long(prev_sd, "SD"), on = c("age_idx", "cohort_idx")
           ][to_long(prev_lb, "LB"), on = c("age_idx", "cohort_idx")
           ][to_long(prev_ub, "UB"), on = c("age_idx", "cohort_idx")]
  dt_prev[, `:=`(age = age_idx - 1, race = race_labels[r], type = "Prevalence")]
  dt_prev <- merge(dt_prev, cohort_lookup, by = "cohort_idx", all.x = TRUE)

  dt_foi <- to_long(foi_mean, "Mean")[to_long(foi_sd, "SD"), on = c("age_idx", "cohort_idx")
          ][to_long(foi_lb, "LB"), on = c("age_idx", "cohort_idx")
          ][to_long(foi_ub, "UB"), on = c("age_idx", "cohort_idx")]
  dt_foi[, `:=`(age = age_idx - 1, race = race_labels[r], type = "Force of Infection")]
  dt_foi <- merge(dt_foi, cohort_lookup, by = "cohort_idx", all.x = TRUE)

  prev_groups[[r]] <- dt_prev
  foi_groups[[r]]  <- dt_foi

  save(dt_prev, dt_foi,
       file = paste0("results/rw_cohort_dt_prev_hat_race_", r, "_",
                      race_short[r], "_RWcohort_expanded.Rdata"))

  cat("Finished predictions for race:", race_labels[r], "\n")
}

df_p_hat_race_final  <- bind_rows(prev_groups)
df_hp_hat_race_final <- bind_rows(foi_groups)


# -----------------------------------------------------------------------------
# 4. Figure: the cohort trajectory itself (rw_alpha / rw_gamma by year)
#    Primary "how did things change across birth cohorts" figure --
#    replaces the old discrete pre/post-1960 comparison. One shared shape
#    (no race facet needed for the trajectory itself), with bin sample
#    counts overlaid so sparsity is visible.
# -----------------------------------------------------------------------------
bin_counts <- as.data.frame(table(cohort_idx)) %>%
  rename(cohort_idx = cohort_idx, n_obs = Freq) %>%
  mutate(cohort_idx = as.integer(as.character(cohort_idx)))

rw_full_summary <- expand_grid(race_idx = race_index, cohort_idx = 1:n_cohorts) %>%
  rowwise() %>%
  mutate(
    alpha_mean = mean(rw_alpha_pop_draws[cohort_idx,] + rw_alpha_dev_draws[race_idx, cohort_idx,]),
    alpha_lo   = quantile(rw_alpha_pop_draws[cohort_idx,] + rw_alpha_dev_draws[race_idx, cohort_idx,], 0.025),
    alpha_hi   = quantile(rw_alpha_pop_draws[cohort_idx,] + rw_alpha_dev_draws[race_idx, cohort_idx,], 0.975),
    gamma_mean = mean(rw_gamma_pop_draws[cohort_idx,] + rw_gamma_dev_draws[race_idx, cohort_idx,]),
    gamma_lo   = quantile(rw_gamma_pop_draws[cohort_idx,] + rw_gamma_dev_draws[race_idx, cohort_idx,], 0.025),
    gamma_hi   = quantile(rw_gamma_pop_draws[cohort_idx,] + rw_gamma_dev_draws[race_idx, cohort_idx,], 0.975)
  ) %>%
  ungroup() %>%
  mutate(race = race_labels[race_idx]) %>%
  left_join(cohort_lookup, by = "cohort_idx")


p_rw_alpha <- ggplot(rw_full_summary, aes(x = birth_cohort, y = alpha_mean, color = race, fill = race)) +
  geom_ribbon(aes(ymin = alpha_lo, ymax = alpha_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Birth cohort year", y = "Cohort effect on alpha (logit scale)",
       title = "Estimated cohort trajectory: max-prevalence (alpha) component, by race") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

p_rw_gamma <- ggplot(rw_full_summary, aes(x = birth_cohort, y = gamma_mean, color = race, fill = race)) +
  geom_ribbon(aes(ymin = gamma_lo, ymax = gamma_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Birth cohort year", y = "Cohort effect on gamma (log scale)",
       title = "Estimated cohort trajectory: rate (gamma) component, by race") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

ggsave("./figures_aug20_2026/rw_alpha_cohort_trajectory.png", p_rw_alpha, width = 10, height = 6)
ggsave("./figures_aug20_2026/rw_gamma_cohort_trajectory.png", p_rw_gamma, width = 10, height = 6)

# -----------------------------------------------------------------------------
# 5. Figure: prevalence-by-age curves, faceted by race, colour = cohort year
#    (replaces the discrete pre/post-1960 two-line comparison, since cohort
#    is now continuous). Plots a representative subset of cohort years for
#    readability.
# -----------------------------------------------------------------------------
cohort_years_to_plot <- cohort_lookup$birth_cohort[
  round(seq(1, n_cohorts, length.out = 10))
]

df_p_curves_subset <- df_p_hat_race_final %>%
  filter(birth_cohort %in% cohort_years_to_plot, age > 0)

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

ggsave("./figures_aug20_2026/predicted_prev_by_cohort_year_RWcohort.png",
       width = 12, height = 8)

df_p_curves_subset %>%
  filter(age <= 10) %>%
  ggplot(aes(x = age, y = Mean, color = birth_cohort, group = birth_cohort)) +
  geom_line(linewidth = 1) +
  facet_wrap(~race, ncol = 3) +
  scale_color_viridis_c(name = "Birth cohort\nyear") +
  xlab("Age") + ylab("Predicted seroprevalence") +
  theme_bw(base_size = 13)
# -----------------------------------------------------------------------------
# 6. Figure: FOI-by-age curves, same structure as above
# -----------------------------------------------------------------------------
df_foi_curves_subset <- df_hp_hat_race_final %>%
  filter(birth_cohort %in% cohort_years_to_plot, age > 0, Mean < 1)

ggplot(df_foi_curves_subset,
       aes(x = age, y = Mean, color = birth_cohort, group = birth_cohort)) +
  geom_line(linewidth = 1) +
  facet_wrap(~race, ncol = 3) +
  scale_color_viridis_c(name = "Birth cohort\nyear") +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  xlab("Age") + ylab("Predicted force of infection") +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "bold"))

ggsave("./figures_aug20_2026/predicted_foi_by_cohort_year_RWcohort_expanded.png",
       width = 12, height = 8)

# -----------------------------------------------------------------------------
# 7. Figure: max prevalence (age 90) by race, at a few reference cohort years
#    Replaces the old pre/post-1960 max-prevalence point-range comparison.
# -----------------------------------------------------------------------------
max_prev_by_race_cohort <- df_p_hat_race_final %>%
  filter(age == fin_age, birth_cohort %in% cohort_years_to_plot)

ggplot(max_prev_by_race_cohort,
       aes(x = race, y = Mean, color = factor(birth_cohort),
           ymin = LB, ymax = UB, group = factor(birth_cohort))) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.7) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  labs(x = "", y = paste0("Max seroprevalence (age ", fin_age, ")"),
       color = "Birth cohort\nyear") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") #bottom

ggsave("./figures_aug20_2026/max_prev_by_race_cohort_RWcohort_expanded.png",
       width = 11, height = 8)

make_max_prev_by_race_cohort <- function(ref_age) {
  df_p_hat_race_final %>%
    filter(age == ref_age, birth_cohort %in% cohort_years_to_plot) %>%
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

# Preserve the original object name for the fin_age (90) case, in case
# anything downstream still refers to `max_prev_by_race_cohort` directly.
max_prev_by_race_cohort <- max_prev_by_race_cohort_list[[paste0("age_", fin_age)]]

for (ref_age in ref_ages) {
  df_this <- max_prev_by_race_cohort_list[[paste0("age_", ref_age)]]
  p_this  <- plot_max_prev_by_race_cohort(df_this, ref_age)

  # Keep the original filename for fin_age (90) so existing downstream
  # references / reports aren't broken; use an age-suffixed name for the
  # new age-5 / age-10 figures.
  out_file <- if (ref_age == fin_age) {
    "./figures_aug20_2026/max_prev_by_race_cohort_RWcohort_expanded.png"
  } else {
    paste0("./figures_aug20_2026/max_prev_by_race_cohort_age", ref_age, "_RWcohort_expanded.png")
  }

  ggsave(out_file, p_this, width = 11, height = 8)
}

# Combined long table across all three reference ages (5, 10, 90) -- one
# save covers every age instead of one save per age.
max_prev_by_race_cohort_all_ages <- bind_rows(max_prev_by_race_cohort_list)
save(max_prev_by_race_cohort_all_ages,
     file = "results/rw_cohort_expanded_max_prev_age5_10_90.Rdata")


# =============================================================================
# 9. Median age of infection (race x cohort_bin), RW-expanded model
#
## =============================================================================

aoi_groups <- list()

for (r in race_index) {

  aoi_mat <- matrix(NA, nrow = n_iter, ncol = n_cohorts)

  for (t in 1:n_cohorts) {
    for (iter in 1:n_iter) {

      g <- exp(par_gamma0_draws[r, iter] + rw_gamma_pop_draws[t, iter] +
                rw_gamma_dev_draws[r, t, iter])

      aoi_mat[iter, t] <- log(2) / g
    }
  }

  dt_aoi <- data.frame(
    cohort_idx = 1:n_cohorts,
    Mean = apply(aoi_mat, 2, mean),
    LB   = apply(aoi_mat, 2, quantile, 0.025),
    UB   = apply(aoi_mat, 2, quantile, 0.975)
  ) %>%
    left_join(cohort_lookup, by = "cohort_idx") %>%
    mutate(race = race_labels[r], type = "Median age of Infection")

  aoi_groups[[r]] <- dt_aoi

  save(dt_aoi,
       file = paste0("results/rw_cohort_dt_aoi_race_", r, "_",
                      race_short[r], "_RWcohort_expanded.Rdata"))

  rm(aoi_mat, dt_aoi)
  cat("Finished median age of infection for race:", race_labels[r], "\n")
}

df_aoi_final <- bind_rows(aoi_groups)

save(df_aoi_final, file = "results/rw_cohort_expanded_df_aoi_all.Rdata")


# -----------------------------------------------------------------------------
# Figure: median age of infection at a few reference cohort years, by race
#    (same reference-year subset used in the prevalence/FOI/max-prev figures,
#    for consistency across the figure set)
# -----------------------------------------------------------------------------
aoi_by_race_cohort <- df_aoi_final %>%
  filter(Mean<100)%>%
  filter(birth_cohort %in% cohort_years_to_plot)

ggplot(aoi_by_race_cohort,
       aes(x = race, y = Mean, color = factor(birth_cohort),
           ymin = LB, ymax = UB, group = factor(birth_cohort))) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.7) +
  labs(x = "", y = "Median age of infection (years)",
       color = "Birth cohort\nyear") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("./figures_aug20_2026/median_age_infection_by_race_cohort_RWcohort_expanded.png",
       width = 11, height = 8)

# -----------------------------------------------------------------------------
# 8. Table/summary: sd_rw_alpha, sd_rw_gamma, theta posteriors
#    (worth reporting directly alongside the trajectory figures)
# -----------------------------------------------------------------------------


param_summary <- tibble(
  parameter = c("sd_rw_alpha_pop", "sd_rw_gamma_pop", "theta"),
  Mean = c(mean(sd_rw_alpha_draws), mean(sd_rw_gamma_draws), mean(theta_draws)),
  LB   = c(quantile(sd_rw_alpha_draws, 0.025), quantile(sd_rw_gamma_draws, 0.025),
           quantile(theta_draws, 0.025)),
  UB   = c(quantile(sd_rw_alpha_draws, 0.975), quantile(sd_rw_gamma_draws, 0.975),
           quantile(theta_draws, 0.975))
)
print(param_summary)
write_csv(param_summary, "results/rw_cohort_key_param_summary.csv")


######

gamma_RW <- exp(par_gamma0_draws[r, ] + rw_gamma_pop_draws[t, ])
gamma_base <- exp(m_gamma0[, r])  # from the base_continuous_bc extraction

cat("RW gamma:   median =", median(gamma_RW), " 95% CI = [",
    quantile(gamma_RW, 0.025), ",", quantile(gamma_RW, 0.975), "]\n")
cat("Base gamma: median =", median(gamma_base), " 95% CI = [",
    quantile(gamma_base, 0.025), ",", quantile(gamma_base, 0.975), "]\n")

# translate to half-life for direct interpretability
cat("RW implied half-life (years):", log(2)/median(gamma_RW), "\n")
cat("Base implied half-life (years):", log(2)/median(gamma_base), "\n")


# *****************************************************************************
#### Visualizations of model convergence ####
# *****************************************************************************

# generate traceplots 
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

m_alpha0 <- mat_all[, paste0("par_alpha0[", 1:n_races, "]"), drop = FALSE]
m_gamma0 <- mat_all[, paste0("par_gamma0[", 1:n_races, "]"), drop = FALSE]

alpha0_labels <- paste0("alpha intercept ", race_labels)
gamma0_labels <- paste0("gamma intercept ", race_labels)

generate_traceplot_function(m_alpha0, alpha0_labels, ncol_facet = 1)
ggsave("./figures_aug20_2026/alpha0_traceplot_RWcohort_expanded.png", width = 12, height = 8)
generate_traceplot_function(m_gamma0, gamma0_labels, ncol_facet = 1)
ggsave("./figures_aug20_2026/gamma0_traceplot_RWcohort_expanded.png", width = 12, height = 8)

generate_densplot_function(m_alpha0, alpha0_labels, ncol_facet = 1)
ggsave("./figures_aug20_2026/alpha0_densplot_RWcohort_expanded.png", width = 12, height = 8)
generate_densplot_function(m_gamma0, gamma0_labels, ncol_facet = 1)
ggsave("./figures_aug20_2026/gamma0_densplot_RWcohort_expanded.png", width = 12, height = 8)

m_alpha_pop <- mat_all[, paste0("rw_alpha_pop[", ref_cohort_bins, "]"), drop = FALSE]
m_gamma_pop <- mat_all[, paste0("rw_gamma_pop[", ref_cohort_bins, "]"), drop = FALSE]

alpha_pop_labels <- paste0("alpha_pop, cohort=", ref_cohort_years)
gamma_pop_labels <- paste0("gamma_pop, cohort=", ref_cohort_years)

generate_traceplot_function(m_alpha_pop, alpha_pop_labels, ncol_facet = 2)
ggsave("./figures_aug20_2026/alpha_pop_traceplot_RWcohort_expanded.png", width = 14, height = 10)
generate_traceplot_function(m_gamma_pop, gamma_pop_labels, ncol_facet = 2)
ggsave("./figures_aug20_2026/gamma_pop_traceplot_RWcohort_expanded.png", width = 14, height = 10)

generate_densplot_function(m_alpha_pop, alpha_pop_labels, ncol_facet = 2)
ggsave("./figures_aug20_2026/alpha_pop_densplot_RWcohort_expanded.png", width = 14, height = 10)
generate_densplot_function(m_gamma_pop, gamma_pop_labels, ncol_facet = 2)
ggsave("./figures_aug20_2026/gamma_pop_densplot_RWcohort_expanded.png", width = 14, height = 10)

for (r in 1:n_races) {
  m_alpha_dev_r <- mat_all[, paste0("rw_alpha_dev[", r, ",", ref_cohort_bins, "]"), drop = FALSE]
  m_gamma_dev_r <- mat_all[, paste0("rw_gamma_dev[", r, ",", ref_cohort_bins, "]"), drop = FALSE]

  alpha_dev_labels <- paste0(race_labels[r], " alpha_dev, cohort=", ref_cohort_years)
  gamma_dev_labels <- paste0(race_labels[r], " gamma_dev, cohort=", ref_cohort_years)

  race_short <- gsub(" ", "", race_labels[r])

  generate_traceplot_function(m_alpha_dev_r, alpha_dev_labels, ncol_facet = 2)
  ggsave(paste0("./figures_aug20_2026/alpha_dev_traceplot_", race_short, "_RWcohort_expanded.png"),
         width = 14, height = 10)
  generate_traceplot_function(m_gamma_dev_r, gamma_dev_labels, ncol_facet = 2)
  ggsave(paste0("./figures_aug20_2026/gamma_dev_traceplot_", race_short, "_RWcohort_expanded.png"),
         width = 14, height = 10)

  generate_densplot_function(m_alpha_dev_r, alpha_dev_labels, ncol_facet = 2)
  ggsave(paste0("./figures_aug20_2026/alpha_dev_densplot_", race_short, "_RWcohort_expanded.png"),
         width = 14, height = 10)
  generate_densplot_function(m_gamma_dev_r, gamma_dev_labels, ncol_facet = 2)
  ggsave(paste0("./figures_aug20_2026/gamma_dev_densplot_", race_short, "_RWcohort_expanded.png"),
         width = 14, height = 10)
}

#all_varnames <- colnames(mat_all)
#study_vars <- grep("^study_re\\[", all_varnames, value = TRUE)
#m_study_re <- mat_all[, study_vars, drop = FALSE]

#if (exists("study_lookup")) {
  #study_idx <- as.integer(gsub("study_re\\[|\\]", "", study_vars))
  #study_labels <- study_lookup$study_abbrev[match(study_idx, study_lookup$study_id_num)]
#} else {
  #study_labels <- study_vars
#}

#generate_traceplot_function(m_study_re, study_labels, ncol_facet = 5)
#ggsave("./figures_aug20_2026/study_re_traceplot_RWexp.png", width = 16, height = 20)
#generate_densplot_function(m_study_re, study_labels, ncol_facet = 5)
#ggsave("./figures_aug20_2026/study_re_densplot_RWexp.png", width = 16, height = 20)



# generate tables for univariate and multivariate convergence statistics

gelman_conv    <- cluster_fit$gelman_conv

gelman_multivariate_PSRF=as.data.frame(gelman_conv$psrf[,1])
colnames(gelman_multivariate_PSRF)="PSRF"
gelman_multivariate_PSRF$param=c("mu_intercept_alpha","mu_intercept_gamma",
"overdispersion_parameter","sd_rw_alpha_pop",
"sd_rw_gamma_pop","sd_rw_dev_alpha",
"sd_rw_dev_gamma","sd_study")

gelman_multivariate_PSRF_table <- gt(cbind.data.frame(parameter=gelman_multivariate_PSRF$param,
PSRF=gelman_multivariate_PSRF$PSRF))
gelman_multivariate_PSRF_table
gtsave(gelman_multivariate_PSRF_table,"./figures_aug20_2026/multivariate_PSRF_RW_cohort_expanded.png")

gelman_par   <- cluster_fit$gelman_par
gelman_univariate_PSRF=as.data.frame(gelman_par$psrf[,1])
colnames(gelman_univariate_PSRF)="PSRF"
gelman_univariate_PSRF$param=c("log_alpha0_H",
"log_alpha0_AIAN","log_alpha0_NHA",
"log_alpha0_NHB","log_alpha0_NHW",
"log_gamma0_H",
"log_gamma0_AIAN","log_gamma0_NHA",
"log_gamma0_NHB","log_gamma0_NHW")

gelman_univariate_PSRF_table <- gt(cbind.data.frame(parameter=gelman_univariate_PSRF$param,
PSRF=gelman_univariate_PSRF$PSRF))
gelman_univariate_PSRF_table
gtsave(gelman_univariate_PSRF_table,"./figures_aug20_2026/univariate_PSRFF_RW_cohort_expanded.png")



