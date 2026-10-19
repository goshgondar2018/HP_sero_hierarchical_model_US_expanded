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

# analytic sample size
total_pop = sum(data_model_filtered$total)

# observed prevalence data points
## cells are pooled by COUNTS rather than by averaging subgroup prevalences.

df_p_curves_obs <- data_model_filtered %>%
  filter(age > 0) %>%
  mutate(era     = factor(ifelse(birth_cohort_cat == 0, "Pre-1960", "1960+"),
                          levels = c("Pre-1960", "1960+"))) %>%
  mutate(study_id_NHANES_indicator=ifelse(study_id=="NHANES1994"|study_id=="NHANES2000",'NHANES','OtherData'))%>%
  group_by(race_fac, era, age = round(age),study_id_NHANES_indicator) %>%
  summarise(n_cell = sum(total, na.rm = TRUE),
            y_cell = sum(HP,    na.rm = TRUE),
            Mean   = y_cell / n_cell,
            LL     = qbeta(0.025, y_cell + 0.5, n_cell - y_cell + 0.5),
            UL     = qbeta(0.975, y_cell + 0.5, n_cell - y_cell + 0.5),
            .groups = "drop") %>%
  mutate(race=ifelse(race_fac=="Hispanic","Hispanic",
ifelse(race_fac=="Non-Hispanic American Indian","AI/AN",
ifelse(race_fac=="Non-Hispanic Asian","API",
ifelse(race_fac=="Non-Hispanic Black", "NHB",
 "NHW")))))

ilogit <- function(logodds) 1 / (1 + exp(-logodds))

n_races <- 5
fin_age <- 90
race_index  <- 1:5
race_labels <- c("Hispanic", "AI/AN", "API", "NHB", "NHW")
race_short  <- c("H", "AIAN", "NHA", "NHB", "NHW")
race_titles <- c("Hispanic",
                 "Non-Hispanic American Indian & Alaska Native",
                 "Non-Hispanic Asian",
                 "Non-Hispanic Black",
                 "Non-Hispanic White")
era_labels  <- c("Pre-1960", "1960+")

dir.create("figures_sep22_2026", showWarnings = FALSE)
dir.create("figures_sep22_2026/binary_cohort_figs", showWarnings = FALSE, recursive = TRUE)
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
    file = "./results_cluster/fit_base.rds",
    kappa_mode = "none", has_study_re = FALSE,
    tag = "Base"
  ),
  kappa_global = list(
    file = "./results_cluster/fit_kappa_global.rds",
    kappa_mode = "global", has_study_re = FALSE,
    tag = "KappaGlobal"
  ),
  kappa_partial = list(
    file = "./results_cluster/fit_kappa_partial.rds",
    kappa_mode = "partial", has_study_re = FALSE,
    tag = "KappaPartial"
  ),
  study_re = list(
    file = "./results_cluster/fit_study_re.rds",
    kappa_mode = "none", has_study_re = TRUE,
    tag = "StudyRE"
  ),
  kappa_partial_study_re = list(
    file = "./results_cluster/fit_kappa_partial_study_re.rds",
    kappa_mode = "partial", has_study_re = TRUE,
    tag = "KappaPartialStudyRE"
    ),
  kappa_global_study_re = list(
    file = "./results_cluster/fit_kappa_global_study_re.rds",
    kappa_mode = "global", has_study_re = TRUE,
    tag = "KappaGlobalStudyRE"
  )
)

# *****************************************************************************
#### Helper functions ####
# *****************************************************************************

# Generalized trace/density plot functions  

generate_traceplot_function <- function(draws, labels, ncol_facet = 1, n_chains = 3) {
  draws_df <- as.data.frame(draws)
  colnames(draws_df) <- labels
  n_iter_per_chain <- nrow(draws_df) / n_chains
  stopifnot(n_iter_per_chain == round(n_iter_per_chain))
  draws_df$iteration <- rep(seq_len(n_iter_per_chain), n_chains)
  draws_df$chain <- rep(factor(seq_len(n_chains)), each = n_iter_per_chain)

  draws_df_long <- pivot_longer(draws_df, cols = -c(iteration,chain),
                                 names_to = "parameter", values_to = "value")
  
  ggplot(draws_df_long, aes(x = iteration, y = value, color = chain)) +
    geom_line(alpha = 0.7,size=1) +
    facet_wrap(~ parameter, scales = "free_y", ncol = ncol_facet) +
    theme_minimal() +
    labs(x = "Iteration", y = "Estimate") +
    theme(legend.position = "bottom",
          legend.title=element_text(size=15),
          legend.text=element_text(size=15),
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

summarize_max_prev_RD <- function(mat, era_label) {
  mat <- as.matrix(mat)
  stopifnot(nrow(mat) == length(race_labels))
  rownames(mat) <- race_labels       # name rows so we don't rely on position

  # subtract every column (iteration)'s NHW row value (referent group)
  rr_mat <- sweep(mat, 2, mat["NHW", ], FUN = "-")

  data.frame(
    Parameters = rownames(rr_mat),
    Mean = apply(rr_mat, 1, mean),
    LB   = apply(rr_mat, 1, quantile, probs = 0.025),
    UB   = apply(rr_mat, 1, quantile, probs = 0.975)
  ) %>%
    mutate(race = Parameters,
           era  = era_label,
           type = "Risk Difference")
}

# Pearcy & Keppel index of disparity, computed PER DRAW.
#
#   ID = 100 * [ sum_i |p_i - p_ref| / n_comparison ] / p_ref
#

summarize_max_prev_IOD <- function(mat, era_label, referent = "NHW") {
  mat <- as.matrix(mat)
  stopifnot(nrow(mat) == length(race_labels))
  rownames(mat) <- race_labels

  dev_abs <- abs(sweep(mat, 2, mat[referent, ], FUN = "-"))
  n_cmp   <- nrow(mat) - 1L
  iod     <- 100 * (colSums(dev_abs) / n_cmp) / mat[referent, ]

  data.frame(
    era  = era_label,
    type = "Index of Disparity",
    Mean = mean(iod),
    LB   = unname(quantile(iod, 0.025)),
    UB   = unname(quantile(iod, 0.975))
  )
}

# change in maximum prevalence between eras, computed PER DRAW
summarize_prev_change <- function(mat_pre, mat_post) {
  mat_pre  <- as.matrix(mat_pre);  rownames(mat_pre)  <- race_labels
  mat_post <- as.matrix(mat_post); rownames(mat_post) <- race_labels

  pct <- 100 * (mat_post - mat_pre) / mat_pre
  abs_pp <- 100 * (mat_post - mat_pre)

  smry <- function(m, lab) data.frame(
    race = rownames(m), type = lab,
    Mean = apply(m, 1, mean),
    LB   = apply(m, 1, quantile, probs = 0.025),
    UB   = apply(m, 1, quantile, probs = 0.975),
    row.names = NULL)

  rbind(smry(abs_pp, "Absolute change (percentage points)"),
        smry(pct,    "Relative change (%)"))
}

summarize_max_prev_RR <- function(mat, era_label) {
  mat <- as.matrix(mat)
  stopifnot(nrow(mat) == length(race_labels))
  rownames(mat) <- race_labels      

  # divide every column (iteration) by that same iteration's NHW row value (referent group)
  rr_mat <- sweep(mat, 2, mat["NHW", ], FUN = "/")

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

      # force of infection: lambda(a) = (dF/da) / (1 - F(a))   [Appendix Eq. 3, 17]
      foi_den <- pmax(1 - p_raw_mat, 1e-10)
      foi_mat <- foi_num / foi_den

      # (gamma*a)^(kappa-1) is infinite at age 0 when kappa < 1, --> FOI undefined
      foi_mat[ages_vec == 0, ] <- NA

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
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/foi_curves_by_race_era_", tag, ".png"), width = 12, height = 8)

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
             aes(x = age, y = Mean, color = era, shape=study_id_NHANES_indicator),
             size = 1.5, alpha = 0.8)+ #, shape = 16
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
  labs(shape="NHANES indicator") +
  theme_bw(base_size = 14) +
  theme(legend.position  = "right",
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
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/predicted_vs_observed_prev_", tag, ".png"), width = 12, height = 8)

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
    theme_bw(base_size = 13) + theme(legend.position = "bottom", strip.text = element_text(face = "bold",size=20),
  axis.text.x=element_text(size=20),
  axis.text.y=element_text(size=20),
  axis.title.x=element_text(size=20),
  axis.title.y=element_text(size=20))
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/prev_curves_by_race_era_", tag, ".png"), width = 12, height = 8)

  max_prev <- df_p_hat_race_final %>% filter(age == fin_age) %>%
    mutate(era = factor(era, levels = era_labels))

  #ggplot(max_prev, aes(x = race, y = Mean, color = race, ymin = LB, ymax = UB, # color = era
                     #   shape = era, group = era)) + # group = era
    #geom_pointrange(size = 1, position = position_dodge(width = 0.6)) +
    #geom_errorbar(width = 0.5, size = 1, position = position_dodge(width = 0.6)) +
    #scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
    #scale_shape_manual(values = c("Pre-1960" = 16, "1960+" = 17)) +
    #xlab("") + ylab("Max seroprevalence (age 90)") +
    #labs(color = "Birth cohort era", shape = "Birth cohort era") +
    #theme_bw(base_size = 16) +
    #theme(legend.position    = "none",#bottom
        #legend.spacing.x   = unit(0.1, "cm"),
        #legend.box         = "horizontal",
        #legend.title       = element_text(size = 14),
        #legend.text        = element_text(size = 14),
        #legend.key.size    = unit(1, "cm"),
       # strip.background   = element_rect(fill = "transparent", color = "transparent"),
  #axis.text.x=element_text(size=20,hjust=1,vjust=1,angle=45),
  #axis.text.y=element_text(size=20),
  #axis.title.x=element_text(size=20),
  #axis.title.y=element_text(size=20))
  
  ggplot(max_prev, aes(x = race, y = Mean, color = race, ymin = LB, ymax = UB, # color = era
                        shape = era)) + # group = era
    geom_pointrange(size = 1, position = position_dodge(width = 0.6)) +
    facet_wrap(~era)+
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
  axis.title.y=element_text(size=20),
  strip.text.x.top =element_text(size=20))
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/max_prev_both_eras_", tag, ".png"), width = 10, height = 8)

  prev_at_5 <- df_p_hat_race_final %>% filter(age == 5) %>%
    mutate(era = factor(era, levels = era_labels))
    
  ggplot(prev_at_5, aes(x = race, y = Mean, color = race, ymin = LB, ymax = UB, # color = era
                        shape = era)) + # group = era
    geom_pointrange(size = 1, position = position_dodge(width = 0.6)) +
    facet_wrap(~era)+
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
  axis.title.y=element_text(size=20),
  strip.text.x.top =element_text(size=20))
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/prev_at_5_both_eras_", tag, ".png"), width = 10, height = 8)

  ## change in max prevalence between eras, computed per posterior draw
  max_prev_percent_change_summarized <- summarize_prev_change(max_prev_draws_pre,
                                                              max_prev_draws_post)
  write_csv(max_prev_percent_change_summarized, paste0("results/percent_change_max_prev_", tag, ".csv")) #_study_re

  ## index of disparity in prev at cutoff of age 90 (with NHW referent),
  ## computed per posterior draw
  max_prev_IOD <- bind_rows(
    summarize_max_prev_IOD(max_prev_draws_pre,  "Pre-1960"),
    summarize_max_prev_IOD(max_prev_draws_post, "1960+")
  ) %>% mutate(era = factor(era, levels = era_labels))

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
  
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/median_age_inf_", tag, ".png"),width = 10, height = 8) #paste0("./figures_sep22_2026/binary_cohort_figs/prev_at_5_both_eras_", tag, ".png"), width = 10, height = 8)


  
  alpha0_labels <- paste0("alpha intercept ", race_labels)
  alpha1_labels <- paste0("alpha slope ", race_labels)
  gamma0_labels <- paste0("gamma intercept ", race_labels)
  gamma1_labels <- paste0("gamma slope ", race_labels)

  generate_traceplot_function(t(m_alpha0), alpha0_labels)
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/alpha0_traceplot_", tag, ".png"), width = 12, height = 8)
  generate_traceplot_function(t(m_alpha1), alpha1_labels)
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/alpha1_traceplot_", tag, ".png"), width = 12, height = 8)
  generate_traceplot_function(t(m_gamma0), gamma0_labels)
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/gamma0_traceplot_", tag, ".png"), width = 12, height = 8)
  generate_traceplot_function(t(m_gamma1), gamma1_labels)
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/gamma1_traceplot_", tag, ".png"), width = 12, height = 8)

  generate_densplot_function(t(m_alpha0), alpha0_labels)
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/alpha0_densplot_", tag, ".png"), width = 12, height = 8)
  generate_densplot_function(t(m_alpha1), alpha1_labels)
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/alpha1_densplot_", tag, ".png"), width = 12, height = 8)
  generate_densplot_function(t(m_gamma0), gamma0_labels)
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/gamma0_densplot_", tag, ".png"), width = 12, height = 8)
  generate_densplot_function(t(m_gamma1), gamma1_labels)
  ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/gamma1_densplot_", tag, ".png"), width = 12, height = 8)

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
    ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/kappa_traceplot_", tag, ".png"), width = 12, height = 8)
    generate_densplot_function(kappa_draws, kappa_labels)
    ggsave(paste0("./figures_sep22_2026/binary_cohort_figs/kappa_densplot_", tag, ".png"), width = 12, height = 8)
  }

  if (cfg$has_study_re && exists("study_lookup")) {
    study_labels <- study_lookup$study_abbrev
  }

  if (!is.null(cluster_fit$gelman_conv)) {
    gelman_multivariate_PSRF <- as.data.frame(cluster_fit$gelman_conv$psrf[, 1])
    colnames(gelman_multivariate_PSRF) <- "PSRF"
    if(model_name=="base") {
      gelman_multivariate_PSRF$param <- c("mu_alpha0", "mu_alpha1","mu_gamma0","mu_gamma1",
                                        "var_alpha0", "var_alpha1", "var_gamma0", "var_gamma1",
                                        "covar_alpha0_1", "covar_alpha0_gamma0", "covar_gamma0_alpha1","covar_gamma1_alpha0","covar_gamma1_alpha1","covar_gamma0_1",
                                      "overdispersion_parameter")
    }
    else if(model_name=='kappa_global') {
      gelman_multivariate_PSRF$param <- c("mu_alpha0", "mu_alpha1","mu_gamma0","mu_gamma1",
                                        "var_alpha0", "var_alpha1", "var_gamma0", "var_gamma1",
                                        "covar_alpha0_1", "covar_alpha0_gamma0", "covar_gamma0_alpha1","covar_gamma1_alpha0","covar_gamma1_alpha1","covar_gamma0_1",
                                      "overdispersion_parameter","kappa")
    }
    else if(model_name=='kappa_partial') {
      gelman_multivariate_PSRF$param <- c("mu_alpha0", "mu_alpha1","mu_gamma0","mu_gamma1",
                                        "var_alpha0", "var_alpha1", "var_gamma0", "var_gamma1",
                                        "covar_alpha0_1", "covar_alpha0_gamma0", "covar_gamma0_alpha1","covar_gamma1_alpha0","covar_gamma1_alpha1","covar_gamma0_1",
                                      "overdispersion_parameter","mu_log_kappa", "sd_log_kappa")
    }
     else if(model_name=='study_re') {
      gelman_multivariate_PSRF$param <- c("mu_alpha0", "mu_alpha1","mu_gamma0","mu_gamma1",
                                        "var_alpha0", "var_alpha1", "var_gamma0", "var_gamma1",
                                        "covar_alpha0_1", "covar_alpha0_gamma0", "covar_gamma0_alpha1","covar_gamma1_alpha0","covar_gamma1_alpha1","covar_gamma0_1",
                                      "overdispersion_parameter","sd_study")
    }
    else if(model_name=='kappa_global_study_re') {
      gelman_multivariate_PSRF$param <- c("mu_alpha0", "mu_alpha1","mu_gamma0","mu_gamma1",
                                        "var_alpha0", "var_alpha1", "var_gamma0", "var_gamma1",
                                        "covar_alpha0_1", "covar_alpha0_gamma0", "covar_gamma0_alpha1","covar_gamma1_alpha0","covar_gamma1_alpha1","covar_gamma0_1",
                                      "overdispersion_parameter","kappa", "sd_study")
    }
    else {
        gelman_multivariate_PSRF$param <- c("mu_alpha0", "mu_alpha1","mu_gamma0","mu_gamma1",
                                        "var_alpha0", "var_alpha1", "var_gamma0", "var_gamma1",
                                        "covar_alpha0_1", "covar_alpha0_gamma0", "covar_gamma0_alpha1","covar_gamma1_alpha0","covar_gamma1_alpha1","covar_gamma0_1",
                                      "overdispersion_parameter","mu_log_kappa", "sd_log_kappa", "sd_study")
     }
    #rownames(cluster_fit$gelman_conv$psrf)
    save_psrf_table(cbind.data.frame(parameter = gelman_multivariate_PSRF$param,
                                      PSRF = gelman_multivariate_PSRF$PSRF),
                    paste0("./figures_sep22_2026/binary_cohort_figs/multivariate_PSRF_", tag, ".png"))
    }
  if (!is.null(cluster_fit$gelman_par)) {
    gelman_univariate_PSRF <- as.data.frame(cluster_fit$gelman_par$psrf[, 1])
    colnames(gelman_univariate_PSRF) <- "PSRF"
    gelman_univariate_PSRF$param <- c("log_alpha0_H", "log_alpha0_AIAN","log_alpha0_API","log_alpha0_NHB","log_alpha0_NHW",
                                      "log_alpha1_H", "log_alpha1_AIAN","log_alpha1_API","log_alpha1_NHB","log_alpha1_NHW",
                                      "log_gamma0_H", "log_gamma0_AIAN","log_gamma0_API","log_gamma0_NHB","log_gamma0_NHW",
                                      "log_gamma1_H", "log_gamma1_AIAN","log_gamma1_API","log_gamma1_NHB","log_gamma1_NHW")
    #rownames(cluster_fit$gelman_par$psrf)
    save_psrf_table(cbind.data.frame(parameter = gelman_univariate_PSRF$param,
                                      PSRF = gelman_univariate_PSRF$PSRF),
                    paste0("./figures_sep22_2026/binary_cohort_figs/univariate_PSRF_", tag, ".png"))
  }
    
  cat("Finished:", model_name, "\n")
  invisible(list(df_p = df_p_hat_race_final, df_foi = df_hp_hat_race_final,
                 df_aoi = df_aoi_final,df_rr = df_max_prev_RR, df_rd = df_max_prev_RD))
}

results_all <- list()
for (model_name in names(model_configs)) {
  results_all[[model_name]] <- process_model(model_name, model_configs[[model_name]])
}

# *****************************************************************************
# Effective sample size, adjusted for beta-binomial overdispersion (theta),
# computed for every model in model_configs (not just the selected one).
#
# For y ~ BetaBinomial(n, p*theta, (1-p)*theta), Var(y) is inflated relative to
# a plain Binomial(n, p) by a design effect DEFF = (n + theta) / (1 + theta).
# The effective sample size is the plain-binomial n that would give the same
# precision: n_eff = n / DEFF = n*(1 + theta) / (n + theta), computed per row
# and per posterior draw of theta, then summed across rows -- same structure
# as total_pop, but downweighted for between-study/subgroup overdispersion.
# The underlying data (and so total_pop) is the same for every model; only
# theta -- and so eff_n -- varies by model.
# *****************************************************************************
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

eff_n_all <- purrr::map_dfr(names(model_configs), function(model_name) {
  cfg    <- model_configs[[model_name]]
  eff_n  <- compute_effective_n(cfg$file, data_model_filtered$total)
  cat(sprintf(
    "Effective sample size (%s, overdispersion-adjusted): %.0f (95%% CrI: %.0f, %.0f); raw total_pop = %.0f\n",
    cfg$tag, eff_n$eff_n_mean, eff_n$eff_n_LB, eff_n$eff_n_UB, eff_n$total_pop
  ))
  cbind(data.frame(model = model_name, tag = cfg$tag), eff_n)
})
write_csv(eff_n_all, "results/effective_sample_size_ALL.csv")



