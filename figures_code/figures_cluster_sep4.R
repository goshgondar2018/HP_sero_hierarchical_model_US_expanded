# *****************************************************************************
#### Load necessary packages ####
# *****************************************************************************

library(tidyverse)
library(data.table)
library(coda)
library(matrixStats)
library(gt)

# *****************************************************************************
#### Data prep  ####
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

df_p_curves_obs <- data_model_filtered %>%
  #filter(age > 0) %>%
  mutate(era     = factor(ifelse(birth_cohort_cat == 0, "Pre-1960", "1960+"),
                          levels = c("Pre-1960", "1960+"))) %>%
  group_by(race_fac, era, age = round(age)) %>%
  summarise(Mean = mean(prev/100, na.rm = TRUE), LL= mean(prev_ci_lb/100, na.rm=TRUE), UL=mean(prev_ci_ub/100, na.rm=TRUE),
            .groups = "drop") %>%
  mutate(race=ifelse(race_fac=="Hispanic","Hispanic",
ifelse(race_fac=="Non-Hispanic American Indian","NH American Indian",
ifelse(race_fac=="Non-Hispanic Asian","NH Asian",
ifelse(race_fac=="Non-Hispanic Black", "NH Black",
 "NH White")))))

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
  study_re = list(
    file = "./cluster_output/fit_study_re.rds",
    kappa_mode = "none", has_study_re = TRUE,
    tag = "StudyRE"
  ),
  kappa_partial_study_re = list(
    file = "./cluster_output/fit_kappa_partial_study_re.rds",
    kappa_mode = "partial", has_study_re = TRUE,
    tag = "KappaPartialStudyRE"
  )
)

# *****************************************************************************
#### Helper functions ####
# *****************************************************************************

# Generalized trace/density plot functions  

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


summarize_max_prev_RD <- function(mat, era_label) {
  mat <- as.matrix(mat)
  stopifnot(nrow(mat) == length(race_labels))
  rownames(mat) <- race_labels       # name rows so we don't rely on position

  # divide every column (iteration) by that same iteration's NH White row value
  rr_mat <- sweep(mat, 2, mat["NH White", ], FUN = "-")

  data.frame(
    Parameters = rownames(rr_mat),
    Mean = apply(rr_mat, 1, mean),
    LB   = apply(rr_mat, 1, quantile, probs = 0.025),
    UB   = apply(rr_mat, 1, quantile, probs = 0.975)
  ) %>%
    mutate(race = Parameters,
           era  = era_label,
           type = "Risk Ratio")
}

summarize_max_prev_RR <- function(mat, era_label) {
  mat <- as.matrix(mat)
  stopifnot(nrow(mat) == length(race_labels))
  rownames(mat) <- race_labels       # name rows so we don't rely on position

  # divide every column (iteration) by that same iteration's NH White row value
  rr_mat <- sweep(mat, 2, mat["NH White", ], FUN = "/")

  data.frame(
    Parameters = rownames(rr_mat),
    Mean = apply(rr_mat, 1, mean),
    LB   = apply(rr_mat, 1, quantile, probs = 0.025),
    UB   = apply(rr_mat, 1, quantile, probs = 0.975)
  ) %>%
    mutate(race = Parameters,
           era  = era_label,
           type = "Risk Ratio")
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

  max_prev_draws_pre  <- matrix(NA_real_, nrow = length(race_index), ncol = n_iter)
  max_prev_draws_post <- matrix(NA_real_, nrow = length(race_index), ncol = n_iter)

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
      #foi_mat[1, ] <- NA

      if (b_idx == 1) {
        max_prev_draws_pre[r, ]  <- p_raw_mat[n_age, ]
      } else {
        max_prev_draws_post[r, ] <- p_raw_mat[n_age, ]
      }

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

  df_max_prev_RR <- bind_rows(
    summarize_max_prev_RR(max_prev_draws_pre,  "Pre-1960"),
    summarize_max_prev_RR(max_prev_draws_post, "1960+")
  ) %>%
    mutate(era = factor(era, levels = era_labels))

  write_csv(df_max_prev_RR, paste0("results/max_prev_RR_", tag, ".csv"))

  df_max_prev_RD <- bind_rows(
    summarize_max_prev_RD(max_prev_draws_pre,  "Pre-1960"),
    summarize_max_prev_RD(max_prev_draws_post, "1960+")
  ) %>%
    mutate(era = factor(era, levels = era_labels))

  write_csv(df_max_prev_RD, paste0("results/max_prev_RD_", tag, ".csv"))

  df_foi_plot <- df_hp_hat_race_final %>%  #filter(age > 0, Mean < 1) %>%
    mutate(era = factor(era, levels = era_labels))

  ggplot(df_foi_plot, aes(x = age, y = Mean, color = era, fill = era, group = era)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2, color = NA) +
    facet_wrap(~race, ncol = 3) +
    scale_color_manual(values = c("Pre-1960" = "#336689", "1960+" = "#E65C4E")) +
    scale_fill_manual(values  = c("Pre-1960" = "#336689", "1960+" = "#E65C4E")) +
    scale_x_continuous(breaks = seq(0, 90, by = 20)) +
    xlab("Age") + ylab("Predicted force of infection") +
    labs(color = "Birth cohort era", fill = "Birth cohort era") +
    theme_bw(base_size = 13) + theme(legend.position = "bottom", strip.text = element_text(face = "bold",size=20),
  axis.legend=element_text(size=30),
  axis.legend.title=element_text(size=30),
  axis.text.x=element_text(size=20),
  axis.text.y=element_text(size=20),
  axis.title.x=element_text(size=20),
  axis.title.y=element_text(size=20))
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/foi_curves_by_race_era_", tag, ".png"), width = 12, height = 8)

# combined predicted and observed prevalence plot
  df_p_plot <- df_p_hat_race_final %>% # filter(age > 0) %>%
    mutate(era = factor(era, levels = era_labels))

  
  ggplot() +
  geom_ribbon(data = df_p_plot,
              aes(x = age, ymin = LB, ymax = UB,
                  fill = era, group = era),
              alpha = 0.2, color = NA) +
  geom_line(data = df_p_plot,
            aes(x = age, y = Mean, color = era, group = era),
            size = 1) +
  geom_point(data = df_p_curves_obs,
             aes(x = age, y = Mean, color = era),
             size = 1, alpha = 0.8, shape = 16)+
  geom_errorbar(data = df_p_curves_obs,
                aes(x = age, ymin = LL, ymax = UL, color = era),
                width = 1.5, alpha = 0.6, size = 0.4) +
  facet_wrap(~race, ncol = 3, scales = "fixed") +
  scale_color_manual(values = c("Pre-1960" = "#336689",
                                "1960+" = "#E65C4E")) +
  scale_fill_manual(values  = c("Pre-1960" = "#336689",
                                "1960+" = "#E65C4E")) +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(0, 1)) +
  xlab("Age") +
  ylab("Seroprevalence") +
  labs(color   = "Birth cohort era",
       fill    = "Birth cohort era") +
  theme_bw(base_size = 14) +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 15),
        legend.title = element_text(size=18),
        strip.background = element_rect(fill  = "transparent",
                                        color = "transparent"),
        strip.text       = element_text(size = 12, face = "bold"),
        axis.text.x      = element_text(size = 15),
        axis.text.y      = element_text(size = 15),
      axis.title.x=element_text(size=15),
      axis.title.y=element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         fill  = guide_legend(override.aes = list(alpha = 0.4)))
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/predicted_vs_observed_prev_", tag, ".png"), width = 12, height = 8)

  #df_p_plot <- df_p_hat_race_final %>% # filter(age > 0) %>%
   # mutate(era = factor(era, levels = era_labels))

  ggplot(df_p_plot, aes(x = age, y = Mean, color = era, fill = era, group = era)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2, color = NA) +
    facet_wrap(~race, ncol = 3) +
    scale_color_manual(values = c("Pre-1960" = "#336689", "1960+" = "#E65C4E")) +
    scale_fill_manual(values  = c("Pre-1960" = "#336689", "1960+" = "#E65C4E")) +
    scale_x_continuous(breaks = seq(0, 90, by = 20)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
    xlab("Age") + ylab("Predicted seroprevalence") +
    labs(color = "Birth cohort era", fill = "Birth cohort era") +
    theme_bw(base_size = 13) + theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
  axis.text.x=element_text(size=20),
  axis.text.y=element_text(size=20),
  axis.title.x=element_text(size=20),
  axis.title.y=element_text(size=20))
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/prev_curves_by_race_era_", tag, ".png"), width = 12, height = 8)

  max_prev <- df_p_hat_race_final %>% filter(age == fin_age) %>%
    mutate(era = factor(era, levels = era_labels))

  ggplot(max_prev, aes(x = race, y = Mean, color = race, ymin = LB, ymax = UB, # color = era
                        shape = era, group = era)) + # group = era
    geom_pointrange(size = 1, position = position_dodge(width = 0.6)) +
    geom_errorbar(width = 0.5, size = 1, position = position_dodge(width = 0.6)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
    scale_shape_manual(values = c("Pre-1960" = 16, "1960+" = 17)) +
    xlab("") + ylab("Max seroprevalence (age 90)") +
    labs(color = "Birth cohort era", shape = "Birth cohort era") +
    theme_bw(base_size = 16) +
    theme(legend.position    = "none",#bottom
        legend.spacing.x   = unit(0.1, "cm"),
        legend.box         = "horizontal",
        legend.title       = element_text(size = 14),
        legend.text        = element_text(size = 14),
        legend.key.size    = unit(1, "cm"),
        strip.background   = element_rect(fill = "transparent", color = "transparent"),
  axis.text.x=element_text(size=20,hjust=1,vjust=1,angle=45),
  axis.text.y=element_text(size=20),
  axis.title.x=element_text(size=20),
  axis.title.y=element_text(size=20))
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/max_prev_both_eras_", tag, ".png"), width = 10, height = 8)

  prev_at_5 <- df_p_hat_race_final %>% filter(age == 5) %>%
    mutate(era = factor(era, levels = era_labels))

  ggplot(prev_at_5, aes(x = race, y = Mean, color = race, ymin = LB, ymax = UB, # color = era
                        shape = era, group = era)) + # group = era
    geom_pointrange(size = 1, position = position_dodge(width = 0.6)) +
    geom_errorbar(width = 0.5, size = 1, position = position_dodge(width = 0.6)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
    scale_shape_manual(values = c("Pre-1960" = 16, "1960+" = 17)) +
    xlab("") + ylab("Max seroprevalence (age 5)") +
    labs(color = "Birth cohort era", shape = "Birth cohort era") +
    theme_bw(base_size = 16) +
    theme(legend.position    = "none",#bottom
        legend.spacing.x   = unit(0.1, "cm"),
        legend.box         = "horizontal",
        legend.title       = element_text(size = 14),
        legend.text        = element_text(size = 14),
        legend.key.size    = unit(1, "cm"),
        strip.background   = element_rect(fill = "transparent", color = "transparent"),
  axis.text.x=element_text(size=20,hjust=1,vjust=1,angle=45),
  axis.text.y=element_text(size=20),
  axis.title.x=element_text(size=20),
  axis.title.y=element_text(size=20))
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/prev_at_5_both_eras_", tag, ".png"), width = 10, height = 8)

  max_prev_percent_change_summarized=max_prev%>%
    group_by(race)%>%
    mutate(mean_percent_change=((Mean[2]-Mean[1])/Mean[1])*100,
         LB_percent_change=((LB[2]-LB[1])/LB[1])*100,
         UB_percent_change=((UB[2]-UB[1])/UB[1])*100)
  write_csv(max_prev_percent_change_summarized, paste0("results/percent_change_max_prev_", tag, ".csv")) #_study_re

  ## index of disparity in prev at cutoff of age 90 (with NHW referent)
  max_prev_IOD=max_prev%>%
     group_by(era)%>%
     mutate(mean_IOD_num=Mean-Mean[5],
         LB_IOD_num=LB-LB[5],
         UB_IOD_num=UB-UB[5])%>%
     group_by(era)%>%
     mutate(mean_IOD_full_num=(sum(mean_IOD_num)/5)*100,
         LB_IOD_full_num=(sum(LB_IOD_num)/5)*100,
        UB_IOD_full_num=(sum(UB_IOD_num)/5)*100)%>%
     group_by(era)%>%
     mutate(mean_IOD=mean_IOD_full_num/Mean[5],
         LB_IOD=mean_IOD_full_num/UB[5],
         UB_IOD=mean_IOD_full_num/LB[5])

  write_csv(max_prev_IOD, file = paste0("results/max_prev_IOD_", tag, ".csv"))
  
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

  df_aoi_final_plot=df_aoi_final%>%
    mutate(era=factor(era, levels = era_labels))

  ggplot(df_aoi_final_plot,
       aes(x = reorder(race,Mean,decreasing=TRUE), y = Mean, 
           ymin = LB, ymax = UB, group = era, color = era, shape = era)) +
    #geom_point()+
    geom_pointrange(size = 1, position = position_dodge(width = 0.6))+
    geom_errorbar(width = 0.5, size = 1, position = position_dodge(width = 0.6)) +
    #facet_wrap(~era, scales = "fixed")+
    theme_bw(base_size = 14) +
    theme(legend.position  = "bottom", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=20),
  legend.text=element_text(size=20),legend.title=element_text(size=20),axis.text.y=element_text(size=20),
        axis.title.x=element_blank(), axis.title.y=element_text(size=20))+ylab("Median age of infection (years)")
  
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/median_age_inf_", tag, ".png"),width = 10, height = 8) #paste0("./figures_aug20_2026/binary_cohort_figs/prev_at_5_both_eras_", tag, ".png"), width = 10, height = 8)


  
  alpha0_labels <- paste0("alpha intercept ", race_labels)
  alpha1_labels <- paste0("alpha slope ", race_labels)
  gamma0_labels <- paste0("gamma intercept ", race_labels)
  gamma1_labels <- paste0("gamma slope ", race_labels)

  generate_traceplot_function(t(m_alpha0), alpha0_labels)
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/alpha0_traceplot_", tag, ".png"), width = 12, height = 8)
  generate_traceplot_function(t(m_alpha1), alpha1_labels)
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/alpha1_traceplot_", tag, ".png"), width = 12, height = 8)
  generate_traceplot_function(t(m_gamma0), gamma0_labels)
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/gamma0_traceplot_", tag, ".png"), width = 12, height = 8)
  generate_traceplot_function(t(m_gamma1), gamma1_labels)
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/gamma1_traceplot_", tag, ".png"), width = 12, height = 8)

  generate_densplot_function(t(m_alpha0), alpha0_labels)
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/alpha0_densplot_", tag, ".png"), width = 12, height = 8)
  generate_densplot_function(t(m_alpha1), alpha1_labels)
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/alpha1_densplot_", tag, ".png"), width = 12, height = 8)
  generate_densplot_function(t(m_gamma0), gamma0_labels)
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/gamma0_densplot_", tag, ".png"), width = 12, height = 8)
  generate_densplot_function(t(m_gamma1), gamma1_labels)
  ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/gamma1_densplot_", tag, ".png"), width = 12, height = 8)

  if (!is.null(m_kappa)) {
    if (cfg$kappa_mode == "global") {
      kappa_draws  <- matrix(m_kappa[1, ], ncol = 1)  # single shared value, all races identical
      kappa_labels <- "kappa (global)"
      kappa_mean=mean(kappa_draws)
      kappa_median=median(kappa_draws)
      kappa_LQ=quantile(kappa_draws,0.025)
      kappa_UQ=quantile(kappa_draws,0.975)
      prob_kappa_greater_than_1=mean(kappa_draws>1)
    } else {
      kappa_draws  <- t(m_kappa)
      kappa_labels <- paste0("kappa ", race_labels)
      kappa_mean=colMeans(kappa_draws)
      kappa_median=colMedians(kappa_draws)
      kappa_LQ=apply(kappa_draws,2,quantile,probs=0.025)
      kappa_UQ=apply(kappa_draws,2,quantile,probs=0.975)
      prob_kappa_greater_than_1=colMeans(kappa_draws>1)
    }
    print(kappa_mean)
    print(as.numeric(kappa_LQ))
    print(as.numeric(kappa_UQ))
    write_csv(cbind.data.frame(mean=kappa_mean,median=kappa_median,LQ=as.numeric(kappa_LQ),UQ=as.numeric(kappa_UQ),
                              prob_greater_than_1=prob_kappa_greater_than_1),
                              paste0("./results/kappa_summary_",tag, ".csv"))

    generate_traceplot_function(kappa_draws, kappa_labels)
    ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/kappa_traceplot_", tag, ".png"), width = 12, height = 8)
    generate_densplot_function(kappa_draws, kappa_labels)
    ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/kappa_densplot_", tag, ".png"), width = 12, height = 8)
  }

  if (cfg$has_study_re && exists("study_lookup")) {
    study_labels <- study_lookup$study_abbrev
    generate_traceplot_function(t(study_re_draws), study_labels, ncol_facet = 5)
    ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/study_re_traceplot_", tag, ".png"), width = 16, height = 20)
    generate_densplot_function(t(study_re_draws), study_labels, ncol_facet = 5)
    ggsave(paste0("./figures_aug20_2026/binary_cohort_figs/study_re_densplot_", tag, ".png"), width = 16, height = 20)
  }

  if (!is.null(cluster_fit$gelman_conv)) {
    gelman_multivariate_PSRF <- as.data.frame(cluster_fit$gelman_conv$psrf[, 1])
    colnames(gelman_multivariate_PSRF) <- "PSRF"
    gelman_multivariate_PSRF$param <- rownames(cluster_fit$gelman_conv$psrf)
    gt_mv <- gt(cbind.data.frame(parameter = gelman_multivariate_PSRF$param,
                                  PSRF = gelman_multivariate_PSRF$PSRF))
    gtsave(gt_mv, paste0("./figures_aug20_2026/binary_cohort_figs/multivariate_PSRF_", tag, ".png"))
  }
  if (!is.null(cluster_fit$gelman_par)) {
    gelman_univariate_PSRF <- as.data.frame(cluster_fit$gelman_par$psrf[, 1])
    colnames(gelman_univariate_PSRF) <- "PSRF"
    gelman_univariate_PSRF$param <- rownames(cluster_fit$gelman_par$psrf)
    gt_uv <- gt(cbind.data.frame(parameter = gelman_univariate_PSRF$param,
                                  PSRF = gelman_univariate_PSRF$PSRF))
    gtsave(gt_uv, paste0("./figures_aug20_2026/binary_cohort_figs/univariate_PSRF_", tag, ".png"))
  }

  cat("Finished:", model_name, "\n")
  invisible(list(df_p = df_p_hat_race_final, df_foi = df_hp_hat_race_final,
                 df_aoi = df_aoi_final,df_rr = df_max_prev_RR, df_rd = df_max_prev_RD))
}

results_all <- list()
for (model_name in names(model_configs)) {
  results_all[[model_name]] <- process_model(model_name, model_configs[[model_name]])
}



