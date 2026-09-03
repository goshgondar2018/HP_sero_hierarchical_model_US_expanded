# *****************************************************************************
#### Load necessary packages ####
# *****************************************************************************

library(tidyverse)
library(data.table)
library(coda)
library(matrixStats)
library(gt)

# *****************************************************************************
#### Data prep (shared across all 4 model variants) ####
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
  #filter(study_id != "SawJRhe2004") %>%
  filter(!(study_id == 'GooJImm2003' & subgroup_n == 384.00)) %>%
  filter(!(study_id == 'FelBirt2007' & subgroup_n == 147.00)) %>%
  filter(!(study_id == 'ParClin2000' & subgroup_n == 2080.00)) %>%
  filter(!(study_id == 'ParClin2000' & subgroup_n == 1040.00)) %>%
  filter(!(study_id == 'CorCanc1990' & age == 0.0410)) %>%
  filter(!(study_id == 'CorCanc1990' & age == 0.2500)) %>%
  filter(!(study_id == 'CorCanc1990' & age == 0.7085))

ilogit <- function(logodds) 1 / (1 + exp(-logodds))

n_races <- 5
fin_age <- 90
race_index  <- 1:5
race_labels <- c("Hispanic", "NH American Indian", "NH Asian", "NH Black", "NH White")
race_short  <- c("H", "AIAN", "NHA", "NHB", "NHW")
era_labels  <- c("Pre-1960", "Post-1960")

dir.create("figures_aug20_2026", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

get_vec <- function(mat_all, base_name, idx) {
  mat_all[, paste0(base_name, "[", idx, "]")]
}
get_par_matrix <- function(mat_all, param_col, n_races) {
  mat_all[, paste0("par[", 1:n_races, ",", param_col, "]"), drop = FALSE]
}

# *****************************************************************************
# Registry of all models
# *****************************************************************************
model_configs <- list(
  base = list(
    file = "./cluster_output/fit_base.rds",
    kappa_mode = "none", has_study_re = FALSE,
    tag = "Base"
  ),
  kappa_global = list(
    file = "./cluster_output/fit_kappa_global.rds",
    kappa_mode = "global", has_study_re = FALSE,
    tag = "KappaGlobal"
  ),
  kappa_partial = list(
    file = "./cluster_output/fit_kappa_partial.rds",
    kappa_mode = "partial", has_study_re = FALSE,
    tag = "KappaPartial"
  ),
  kappa_partial_study_re = list(
    file = "./cluster_output/fit_kappa_partial_study_re.rds",
    kappa_mode = "partial", has_study_re = TRUE,
    tag = "KappaPartialStudyRE"
  )
)

# *****************************************************************************
#### Generalized trace/density plot functions  ####
# *****************************************************************************
generate_traceplot_function <- function(draws, labels, ncol_facet = 1) {
  draws_df <- as.data.frame(draws)
  colnames(draws_df) <- labels
  draws_df$iteration <- seq_len(nrow(draws_df))
  draws_df_long <- pivot_longer(draws_df, cols = -iteration,
                                 names_to = "parameter", values_to = "value")
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
  draws_df_long <- pivot_longer(draws_df, cols = -iteration,
                                 names_to = "parameter", values_to = "value")
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

# *****************************************************************************
#### Process one model: prediction curves, max prev, AoI, trace/dens, PSRF ####
# *****************************************************************************
process_model <- function(model_name, cfg) {

  cat("\n===== Processing:", model_name, "=====\n")

  cluster_fit <- readRDS(cfg$file)
  mat_all     <- as.matrix(cluster_fit$mcmc_all)
  n_iter      <- nrow(mat_all)
  tag         <- cfg$tag

  m_alpha0 <- t(get_par_matrix(mat_all, 1, n_races))  
  m_alpha1 <- t(get_par_matrix(mat_all, 2, n_races))
  m_gamma0 <- t(get_par_matrix(mat_all, 3, n_races))
  m_gamma1 <- t(get_par_matrix(mat_all, 4, n_races))

  m_kappa <- NULL
  if (cfg$kappa_mode == "global") {
    kappa_vec_shared <- mat_all[, "kappa"]                       
    m_kappa <- matrix(kappa_vec_shared, nrow = n_races, ncol = n_iter, byrow = TRUE)
  } else if (cfg$kappa_mode == "partial") {
    m_kappa <- t(mat_all[, paste0("kappa[", 1:n_races, "]"), drop = FALSE])  
  }

  study_re_draws <- NULL
  if (cfg$has_study_re) {
    n_studies <- length(unique(as.numeric(as.factor(data_model_filtered$study_id))))
    study_re_draws <- t(sapply(1:n_studies, function(s) get_vec(mat_all, "study_re", s)))
  }
  
  theta_draws <- mat_all[, "theta"]

  ages_vec <- 0:fin_age
  n_age    <- fin_age + 1

  prev_groups <- vector("list", length(race_index))
  foi_groups  <- vector("list", length(race_index))

  for (r in race_index) {

    prev_mean <- prev_sd <- prev_lb <- prev_ub <- matrix(NA_real_, n_age, 2)
    foi_mean  <- foi_sd  <- foi_lb  <- foi_ub  <- matrix(NA_real_, n_age, 2)

    for (b_idx in 1:2) {
      b <- b_idx - 1  

      a_vec <- ilogit(m_alpha0[r, ] + m_alpha1[r, ] * b)
      g_vec <- exp(   m_gamma0[r, ] + m_gamma1[r, ] * b)

      g_age_raw <- outer(ages_vec, g_vec, FUN = "*")  

      a_rep <- rep(a_vec, each = n_age)
      g_rep <- rep(g_vec, each = n_age)

      if (cfg$kappa_mode == "none") {
        exp_neg       <- exp(-g_age_raw)
        one_minus_exp <- 1 - exp_neg
        p_raw_mat     <- one_minus_exp * a_rep
        p_mat         <- 0.001 + 0.998 * p_raw_mat
        foi_num       <- a_rep * g_rep * exp_neg

      } else {
        kappa_vec     <- m_kappa[r, ]
        kappa_rep_mat <- matrix(kappa_vec, nrow = n_age, ncol = length(g_vec), byrow = TRUE)
        g_age_k       <- g_age_raw^kappa_rep_mat
        exp_neg       <- exp(-g_age_k)
        one_minus_exp <- 1 - exp_neg
        p_raw_mat     <- one_minus_exp * a_rep
        p_mat         <- 0.001 + 0.998 * p_raw_mat

        kappa_rep <- rep(kappa_vec, each = n_age)
        g_age_raw_vec <- as.vector(g_age_raw)
        foi_num <- a_rep * g_rep * kappa_rep * (g_age_raw_vec)^(kappa_rep - 1) * exp_neg
      }

      foi_den <- 1 - p_mat * one_minus_exp + 1e-10
      foi_mat <- foi_num / foi_den
      foi_mat[1, ] <- NA

      prev_mean[, b_idx] <- rowMeans2(p_raw_mat, na.rm = TRUE)
      prev_sd[, b_idx]   <- rowSds(p_raw_mat,   na.rm = TRUE)
      pq <- rowQuantiles(p_raw_mat, probs = c(0.025, 0.975), na.rm = TRUE)
      prev_lb[, b_idx] <- pq[, 1]; prev_ub[, b_idx] <- pq[, 2]

      foi_mean[, b_idx] <- rowMeans2(foi_mat, na.rm = TRUE)
      foi_sd[, b_idx]   <- rowSds(foi_mat,   na.rm = TRUE)
      fq <- rowQuantiles(foi_mat, probs = c(0.025, 0.975), na.rm = TRUE)
      foi_lb[, b_idx] <- fq[, 1]; foi_ub[, b_idx] <- fq[, 2]
    }

    to_long <- function(mat, stat_name) {
      as.data.table(reshape2::melt(mat, varnames = c("age_idx", "era_idx"),
                                    value.name = stat_name))
    }

    dt_prev <- to_long(prev_mean, "Mean")[to_long(prev_sd, "SD"), on = c("age_idx", "era_idx")
             ][to_long(prev_lb, "LB"), on = c("age_idx", "era_idx")
             ][to_long(prev_ub, "UB"), on = c("age_idx", "era_idx")]
    dt_prev[, `:=`(age = age_idx - 1, race = race_labels[r], type = "Prevalence",
                    era = era_labels[era_idx])]

    dt_foi <- to_long(foi_mean, "Mean")[to_long(foi_sd, "SD"), on = c("age_idx", "era_idx")
            ][to_long(foi_lb, "LB"), on = c("age_idx", "era_idx")
            ][to_long(foi_ub, "UB"), on = c("age_idx", "era_idx")]
    dt_foi[, `:=`(age = age_idx - 1, race = race_labels[r], type = "Force of Infection",
                   era = era_labels[era_idx])]

    prev_groups[[r]] <- dt_prev
    foi_groups[[r]]  <- dt_foi
  }

  df_p_hat_race_final  <- bind_rows(prev_groups)
  df_hp_hat_race_final <- bind_rows(foi_groups)

  save(df_p_hat_race_final, df_hp_hat_race_final,
       file = paste0("results/", tag, "_prev_foi_all.Rdata"))

  df_foi_plot <- df_hp_hat_race_final %>% filter(age > 0, Mean < 1) %>%
    mutate(era = factor(era, levels = era_labels))

  ggplot(df_foi_plot, aes(x = age, y = Mean, color = era, fill = era, group = era)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2, color = NA) +
    facet_wrap(~race, ncol = 3) +
    scale_color_manual(values = c("Pre-1960" = "#336689", "Post-1960" = "#E65C4E")) +
    scale_fill_manual(values  = c("Pre-1960" = "#336689", "Post-1960" = "#E65C4E")) +
    scale_x_continuous(breaks = seq(0, 90, by = 20)) +
    xlab("Age") + ylab("Predicted force of infection") +
    labs(color = "Birth cohort era", fill = "Birth cohort era") +
    theme_bw(base_size = 13) + theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
  ggsave(paste0("./figures_aug20_2026/foi_curves_by_race_era_", tag, ".png"), width = 12, height = 8)

  df_p_plot <- df_p_hat_race_final %>% filter(age > 0) %>%
    mutate(era = factor(era, levels = era_labels))

  ggplot(df_p_plot, aes(x = age, y = Mean, color = era, fill = era, group = era)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2, color = NA) +
    facet_wrap(~race, ncol = 3) +
    scale_color_manual(values = c("Pre-1960" = "#336689", "Post-1960" = "#E65C4E")) +
    scale_fill_manual(values  = c("Pre-1960" = "#336689", "Post-1960" = "#E65C4E")) +
    scale_x_continuous(breaks = seq(0, 90, by = 20)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
    xlab("Age") + ylab("Predicted seroprevalence") +
    labs(color = "Birth cohort era", fill = "Birth cohort era") +
    theme_bw(base_size = 13) + theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
  ggsave(paste0("./figures_aug20_2026/prev_curves_by_race_era_", tag, ".png"), width = 12, height = 8)

  max_prev <- df_p_hat_race_final %>% filter(age == fin_age) %>%
    mutate(era = factor(era, levels = era_labels))

  ggplot(max_prev, aes(x = race, y = Mean, color = era, ymin = LB, ymax = UB,
                        shape = era, group = era)) +
    geom_pointrange(size = 1, position = position_dodge(width = 0.6)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
    scale_shape_manual(values = c("Pre-1960" = 16, "Post-1960" = 17)) +
    xlab("") + ylab("Max seroprevalence (age 90)") +
    labs(color = "Birth cohort era", shape = "Birth cohort era") +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0("./figures_aug20_2026/max_prev_both_eras_", tag, ".png"), width = 10, height = 8)

  aoi_groups <- list()
  for (r in race_index) {
    aoi_mat <- matrix(NA, nrow = n_iter, ncol = 2)
    for (b_idx in 1:2) {
      b <- b_idx - 1
      g <- exp(m_gamma0[r, ] + m_gamma1[r, ] * b)
      if (cfg$kappa_mode == "none") {
        aoi_mat[, b_idx] <- log(2) / g
      } else {
        kappa_vec <- m_kappa[r, ]
        aoi_mat[, b_idx] <- (log(2))^(1/kappa_vec) / g
      }
    }
    dt_aoi <- data.frame(
      era_idx = 1:2,
      Mean = apply(aoi_mat, 2, mean),
      LB   = apply(aoi_mat, 2, quantile, 0.025),
      UB   = apply(aoi_mat, 2, quantile, 0.975)
    ) %>% mutate(era = era_labels[era_idx], race = race_labels[r],
                  type = "Median age of Infection")
    aoi_groups[[r]] <- dt_aoi
  }
  df_aoi_final <- bind_rows(aoi_groups)
  write_csv(df_aoi_final, paste0("results/median_age_inf_", tag, ".csv"))

  alpha0_labels <- paste0("alpha intercept ", race_labels)
  alpha1_labels <- paste0("alpha slope ", race_labels)
  gamma0_labels <- paste0("gamma intercept ", race_labels)
  gamma1_labels <- paste0("gamma slope ", race_labels)

  generate_traceplot_function(t(m_alpha0), alpha0_labels)
  ggsave(paste0("./figures_aug20_2026/alpha0_traceplot_", tag, ".png"), width = 12, height = 8)
  generate_traceplot_function(t(m_alpha1), alpha1_labels)
  ggsave(paste0("./figures_aug20_2026/alpha1_traceplot_", tag, ".png"), width = 12, height = 8)
  generate_traceplot_function(t(m_gamma0), gamma0_labels)
  ggsave(paste0("./figures_aug20_2026/gamma0_traceplot_", tag, ".png"), width = 12, height = 8)
  generate_traceplot_function(t(m_gamma1), gamma1_labels)
  ggsave(paste0("./figures_aug20_2026/gamma1_traceplot_", tag, ".png"), width = 12, height = 8)

  generate_densplot_function(t(m_alpha0), alpha0_labels)
  ggsave(paste0("./figures_aug20_2026/alpha0_densplot_", tag, ".png"), width = 12, height = 8)
  generate_densplot_function(t(m_alpha1), alpha1_labels)
  ggsave(paste0("./figures_aug20_2026/alpha1_densplot_", tag, ".png"), width = 12, height = 8)
  generate_densplot_function(t(m_gamma0), gamma0_labels)
  ggsave(paste0("./figures_aug20_2026/gamma0_densplot_", tag, ".png"), width = 12, height = 8)
  generate_densplot_function(t(m_gamma1), gamma1_labels)
  ggsave(paste0("./figures_aug20_2026/gamma1_densplot_", tag, ".png"), width = 12, height = 8)

  if (!is.null(m_kappa)) {
    if (cfg$kappa_mode == "global") {
      kappa_draws  <- matrix(m_kappa[1, ], ncol = 1)  # single shared value, all races identical
      kappa_labels <- "kappa (global)"
    } else {
      kappa_draws  <- t(m_kappa)
      kappa_labels <- paste0("kappa ", race_labels)
    }
    generate_traceplot_function(kappa_draws, kappa_labels)
    ggsave(paste0("./figures_aug20_2026/kappa_traceplot_", tag, ".png"), width = 12, height = 8)
    generate_densplot_function(kappa_draws, kappa_labels)
    ggsave(paste0("./figures_aug20_2026/kappa_densplot_", tag, ".png"), width = 12, height = 8)
  }

  if (cfg$has_study_re && exists("study_lookup")) {
    study_labels <- study_lookup$study_abbrev
    generate_traceplot_function(t(study_re_draws), study_labels, ncol_facet = 5)
    ggsave(paste0("./figures_aug20_2026/study_re_traceplot_", tag, ".png"), width = 16, height = 20)
    generate_densplot_function(t(study_re_draws), study_labels, ncol_facet = 5)
    ggsave(paste0("./figures_aug20_2026/study_re_densplot_", tag, ".png"), width = 16, height = 20)
  }

  if (!is.null(cluster_fit$gelman_conv)) {
    gelman_multivariate_PSRF <- as.data.frame(cluster_fit$gelman_conv$psrf[, 1])
    colnames(gelman_multivariate_PSRF) <- "PSRF"
    gelman_multivariate_PSRF$param <- rownames(cluster_fit$gelman_conv$psrf)
    gt_mv <- gt(cbind.data.frame(parameter = gelman_multivariate_PSRF$param,
                                  PSRF = gelman_multivariate_PSRF$PSRF))
    gtsave(gt_mv, paste0("./figures_aug20_2026/multivariate_PSRF_", tag, ".png"))
  }
  if (!is.null(cluster_fit$gelman_par)) {
    gelman_univariate_PSRF <- as.data.frame(cluster_fit$gelman_par$psrf[, 1])
    colnames(gelman_univariate_PSRF) <- "PSRF"
    gelman_univariate_PSRF$param <- rownames(cluster_fit$gelman_par$psrf)
    gt_uv <- gt(cbind.data.frame(parameter = gelman_univariate_PSRF$param,
                                  PSRF = gelman_univariate_PSRF$PSRF))
    gtsave(gt_uv, paste0("./figures_aug20_2026/univariate_PSRF_", tag, ".png"))
  }

  cat("Finished:", model_name, "\n")
  invisible(list(df_p = df_p_hat_race_final, df_foi = df_hp_hat_race_final, df_aoi = df_aoi_final))
}

results_all <- list()
for (model_name in names(model_configs)) {
  results_all[[model_name]] <- process_model(model_name, model_configs[[model_name]])
}

