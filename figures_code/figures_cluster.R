# *****************************************************************************
#### Load necessary packages ####
# *****************************************************************************
library(ggplot2)
library(coda)
library(data.table)

# *****************************************************************************
#### Data prep ####
# *****************************************************************************
 
## input data 
load("data/df_hp_usa_w_final.Rdata")

data_model <- df_hp_usa_w_final %>%
  arrange(race_f, year)%>%
  filter(race_f != 'Other')%>%
  mutate(study_id_fac = factor(study_id),
         race_fac = factor(race_f),
         time = year - 1988,
         bc_scaled = birth_cohort - min(birth_cohort))%>%
  mutate(birth_cohort_cat=ifelse(birth_cohort<1960,0,1))

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

## model output
cluster_fit <- readRDS("./cluster_output/fit_study_re.Rds") #readRDS("./cluster_output/fit_study_re.Rds")
mcmc_all    <- cluster_fit$mcmc_all
mat_all <- as.matrix(mcmc_all)
cat("Available columns (sample):\n")
print(head(colnames(mat_all), 20))

n_races <- 5   
fin_age <- 90
n_iter <- nrow(mat_all)

# *****************************************************************************
#### Helper functions ####
# *****************************************************************************

ilogit <- function(logodds) 1 / (1 + exp(-logodds))
logit  <- function(p) log(p / (1 - p))

## melt 3D array, summarize, rename, recode Var3 to 0/1
summarize_predictions <- function(arr) {
  dt <- as.data.table(reshape2::melt(arr))
  dt <- dt[, .(Mean = mean(value,    na.rm = TRUE),
               SD   = sd(value,      na.rm = TRUE),
               LB   = quantile(value, 0.025, na.rm = TRUE),
               UB   = quantile(value, 0.975, na.rm = TRUE)),
           by = .(Var1, Var3)]
  setnames(dt, c("Var1", "Var3"), c("age", "birth_cohort_cat"))
  dt[, birth_cohort_cat := birth_cohort_cat - 1]  # recode 1/2 → 0/1
  dt
}

## label binary cohort era values
add_era_label <- function(dt) {
  dt %>% mutate(birth_cohort_label = ifelse(birth_cohort_cat == 0, "Pre-1960", "1960+"))
}

## reconstruct m_alpha0, m_alpha1, m_gamma0, m_gamma1 as [n_iter x n_races] matrices 
get_par_matrix <- function(mat_all, param_col, n_races) {
  cols <- paste0("par[", 1:n_races, ",", param_col, "]")
  out <- mat_all[, cols, drop = FALSE]
  colnames(out) <- paste0("param", 1:n_races)  
  out
}

## process results output
m_alpha0 <- get_par_matrix(mat_all, 1, n_races)
m_alpha1 <- get_par_matrix(mat_all, 2, n_races)
m_gamma0 <- get_par_matrix(mat_all, 3, n_races)
m_gamma1 <- get_par_matrix(mat_all, 4, n_races)

m_alpha0_mub <- mat_all[, "mub[1]", drop = FALSE]
m_alpha1_mub <- mat_all[, "mub[2]", drop = FALSE]
m_gamma0_mub <- mat_all[, "mub[3]", drop = FALSE]
m_gamma1_mub <- mat_all[, "mub[4]", drop = FALSE]


# *****************************************************************************
#### Prediction loop — all races ####
# *****************************************************************************

race_index <- 1:5
race_labels <- c("Hispanic", "NH American Indian", "NH Asian", "NH Black", "NH White")
race_short  <- c("H", "AIAN", "NHA", "NHB", "NHW")
race_titles <- c("Hispanic",
                 "Non-Hispanic American Indian & Alaska Native",
                 "Non-Hispanic Asian",
                 "Non-Hispanic Black",
                 "Non-Hispanic White")

prev_groups <- list()
foi_groups  <- list()

for (r in race_index) {
  
  alpha_est <- array(NA, dim = c(n_iter, 2))
  gamma_est <- array(NA, dim = c(n_iter, 2))
  p_hat     <- array(NA, dim = c(fin_age + 1, n_iter, 2))
  hp_hat    <- array(NA, dim = c(fin_age + 1, n_iter, 2))
  
  for (iter in 1:n_iter) {
    for (b in 0:1) {
      a <- ilogit(m_alpha0[iter, r] + m_alpha1[iter, r] * b)
      g <- exp(  m_gamma0[iter, r] + m_gamma1[iter, r] * b)
      
      alpha_est[iter, b + 1] <- a
      gamma_est[iter, b + 1] <- g
      
      for (age in 0:fin_age) {
        p_raw <- a * (1 - exp(-g * age))
        p     <- 0.001 + 0.998 * p_raw
        
        p_hat[age + 1, iter, b + 1] <- p_raw
        
        if(age == 0){
          hp_hat[age + 1, iter, b + 1] <- NA
        } else {
          hp_hat[age + 1, iter, b + 1] <- (a * g * exp(-g * age)) /
                                           (1 - p * (1 - exp(-g * age)) + 1e-10)
        }

      }
    }
  }
  
  dt_prev <- summarize_predictions(p_hat)
  dt_prev[, `:=`(race = race_labels[r], type = "Prevalence")]
  
  dt_foi  <- summarize_predictions(hp_hat)
  dt_foi[,  `:=`(race = race_labels[r], type = "Force of Infection (λ)")]
  
  prev_groups[[r]] <- dt_prev
  foi_groups[[r]]  <- dt_foi
  
  rm(alpha_est, gamma_est, p_hat, hp_hat, dt_prev, dt_foi) 
}


# *****************************************************************************
#### Final joined datasets ####
# *****************************************************************************

df_p_hat_race_final  <- bind_rows(lapply(prev_groups, as.data.frame))
df_hp_hat_race_final <- bind_rows(lapply(foi_groups,  as.data.frame))

save(df_p_hat_race_final, df_hp_hat_race_final,
     file = "results/all_predicted_prev_FOI_data_study_re.Rdata") # _study_re


# *****************************************************************************
#### Max prevalence at multiple ages ####
# *****************************************************************************

v_race <- c("Hispanic", "NH American Indian", "NH Asian", "NH Black", "NH White")

# summarize race-specific estimates for both eras
summarize_max_prev <- function(mat, era_label) {
  as.data.frame(mat) %>%
    gather(key = "Parameters", value = "Value") %>%
    group_by(Parameters) %>%
    summarise(Mean = mean(Value),
              LB   = quantile(Value, 0.025),
              UB   = quantile(Value, 0.975),
              .groups = "drop") %>%
    mutate(race = v_race,
           era  = era_label,
           type = "Max Prevalence")
}

# summarize population mean reference lines for both eras
summarize_mub_mean <- function(mat, era_label) {
  data.frame(max_prev = c(mat)) %>%
    summarise(Mean = mean(max_prev),
              LB   = quantile(max_prev, 0.025),
              UB   = quantile(max_prev, 0.975)) %>%
    mutate(era = era_label)
}

# max-prevalence calculator for any max age (age_cutoff)
compute_max_prev_by_age <- function(age_cutoff) {

  # pre-1960 (bc_binary = 0)
  max_prev_pre     <- ilogit(m_alpha0) * (1 - exp(-exp(m_gamma0) * age_cutoff))
  max_prev_pre_mub <- ilogit(m_alpha0_mub) * (1 - exp(-exp(m_gamma0_mub) * age_cutoff))

  # 1960+ (bc_binary = 1)
  max_prev_post     <- ilogit(m_alpha0 + m_alpha1) * (1 - exp(-exp(m_gamma0 + m_gamma1) * age_cutoff))
  max_prev_post_mub <- ilogit(m_alpha0_mub + m_alpha1_mub) * (1 - exp(-exp(m_gamma0_mub + m_gamma1_mub) * age_cutoff))

  df_max_prev_both <- bind_rows(
    summarize_max_prev(max_prev_pre,  "Pre-1960"),
    summarize_max_prev(max_prev_post, "1960+")
  ) %>%
    mutate(era = factor(era, levels = c("Pre-1960", "1960+")),
           age_cutoff = age_cutoff)

  df_max_prev_mean_both <- bind_rows(
    summarize_mub_mean(max_prev_pre_mub,  "Pre-1960"),
    summarize_mub_mean(max_prev_post_mub, "1960+")
  ) %>%
    reshape2::melt(id.vars = "era", variable.name = "variables", value.name = "FE") %>%
    mutate(Type = ifelse(variables == "Mean", "Overall mean", "Credible Interval"),
           age_cutoff = age_cutoff)

  list(df_max_prev_both = df_max_prev_both,
       df_max_prev_mean_both = df_max_prev_mean_both)
}

# cutoffs at age 5 and age 10
max_prev_ages   <- c(5, 10, fin_age)
max_prev_by_age <- lapply(max_prev_ages, compute_max_prev_by_age)
names(max_prev_by_age) <- paste0("age_", max_prev_ages)
df_max_prev_both      <- max_prev_by_age[[paste0("age_", fin_age)]]$df_max_prev_both
df_max_prev_mean_both <- max_prev_by_age[[paste0("age_", fin_age)]]$df_max_prev_mean_both
df_max_prev_age5  <- max_prev_by_age$age_5$df_max_prev_both
df_max_prev_age10 <- max_prev_by_age$age_10$df_max_prev_both
df_max_prev_all_ages <- bind_rows(lapply(max_prev_by_age, `[[`, "df_max_prev_both"))
save(df_max_prev_all_ages, file = "./results/max_prev_age5_10_90_study_re.Rdata")

## plot for cutoff at age 90
ggplot(df_max_prev_both,
       aes(x = race, y = Mean, color = race,
           ymin = LB, ymax = UB,
           shape = era, group = era)) +
  #geom_point()+
  geom_pointrange(size = 1, position = position_dodge(width = 0.6))+
  geom_errorbar(width = 0.5, size = 1, position = position_dodge(width = 0.6)) +
  #geom_hline(data = df_max_prev_mean_both,
             #aes(yintercept = FE, linetype = interaction(Type, era),
                # color = NULL),
            # col = "#336689", size = 0.8) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  scale_shape_manual(values = c("Pre-1960" = 16, "1960+" = 17)) +
  scale_linetype_manual(
    values = c("Overall mean.Pre-1960"    = "solid",
               "Credible Interval.Pre-1960" = "dashed",
               "Overall mean.1960+"   = "solid",
               "Credible Interval.1960+" = "dashed")) +
  #scale_color_manual(values = colors_5) +
  xlab("") + ylab("Max seroprevalence (age 90)") +
  labs(shape = "Birth cohort era") +
  theme_bw(base_size = 20) +
  theme(legend.position    = "none",#bottom
        legend.spacing.x   = unit(0.1, "cm"),
        legend.box         = "horizontal",
        legend.title       = element_text(size = 14),
        legend.text        = element_text(size = 14),
        legend.key.size    = unit(1, "cm"),
        strip.background   = element_rect(fill = "transparent", color = "transparent"),
        axis.text.x        = element_text(hjust = 1, vjust = 1, angle = 45)) +
  guides(color    = "none",
         shape    = guide_legend(override.aes = list(size = 3)),
         linetype = "none")  # hline legend gets busy — suppress

ggsave("./figures_aug20_2026/max_prev_both_eras_study_re.png", width = 10, height = 8) #_study_re

## plot for cutoff at age 5
ggplot(df_max_prev_age5,
       aes(x = race, y = Mean, color = race,
           ymin = LB, ymax = UB,
           shape = era, group = era)) +
  #geom_point()+
  geom_pointrange(size = 1, position = position_dodge(width = 0.6))+
  geom_errorbar(width = 0.5, size = 1, position = position_dodge(width = 0.6)) +
  #geom_hline(data = df_max_prev_mean_both,
             #aes(yintercept = FE, linetype = interaction(Type, era),
                # color = NULL),
            # col = "#336689", size = 0.8) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  scale_shape_manual(values = c("Pre-1960" = 16, "1960+" = 17)) +
  scale_linetype_manual(
    values = c("Overall mean.Pre-1960"    = "solid",
               "Credible Interval.Pre-1960" = "dashed",
               "Overall mean.1960+"   = "solid",
               "Credible Interval.1960+" = "dashed")) +
  #scale_color_manual(values = colors_5) +
  xlab("") + ylab("Max seroprevalence (age 5)") +
  labs(shape = "Birth cohort era") +
  theme_bw(base_size = 20) +
  theme(legend.position    = "bottom",
        legend.spacing.x   = unit(0.1, "cm"),
        legend.box         = "horizontal",
        legend.title       = element_text(size = 14),
        legend.text        = element_text(size = 14),
        legend.key.size    = unit(1, "cm"),
        strip.background   = element_rect(fill = "transparent", color = "transparent"),
        axis.text.x        = element_text(hjust = 1, vjust = 1, angle = 45)) +
  guides(color    = "none",
         shape    = guide_legend(override.aes = list(size = 3)),
         linetype = "none")  # hline legend gets busy — suppress

ggsave("./figures_aug20_2026/max_prev_age5_both_eras_study_re.png", width = 10, height = 8) #_study_re

## summarize % changes in max prevalence across birth cohorts by race/ethnic group
df_max_prev_both_summarized=df_max_prev_both%>%
  group_by(race)%>%
  mutate(mean_percent_change=((Mean[2]-Mean[1])/Mean[1])*100,
         LB_percent_change=((LB[2]-LB[1])/LB[1])*100,
         UB_percent_change=((UB[2]-UB[1])/UB[1])*100)
save(df_max_prev_both_summarized, file = "./results/percent_change_max_prev_study_re.Rdata") #_study_re

## index of disparity in prev at cutoff of age 90 (with NHW referent)
df_max_prev_IOD=df_max_prev_both%>%
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

save(df_max_prev_IOD, file = "./results/prev_IOD_study_re.Rdata") #_study_re
  
## compute ratios of max prev by age cutoff of 90 for each race/ethnicity compared to NHW referent
df_max_prev_RR=df_max_prev_both%>%
  group_by(era)%>%
  mutate(mean_RR=Mean/Mean[5],
         LB_RR=UB/UB[5],
         UB_RR=LB/LB[5])

save(df_max_prev_RR, file = "./results/prev_RR_comparisons_study_re.Rdata") #_study_re

# *****************************************************************************
#### FOI and prevalence curves across ages by race and era ####
# *****************************************************************************

## FOI plot
df_foi_curves <- df_hp_hat_race_final %>%
  mutate(era = factor(ifelse(birth_cohort_cat == 0, "Pre-1960", "1960+"),
                      levels = c("Pre-1960", "1960+")))

ggplot(df_foi_curves,
       aes(x = age, y = Mean, color = era, fill = era, group = era)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2, color = NA) +
  facet_wrap(~race, ncol = 3) +
  scale_color_manual(values = c("Pre-1960" = "#336689",
                                "1960+" = "#E65C4E")) +
  scale_fill_manual(values  = c("Pre-1960" = "#336689",
                                "1960+" = "#E65C4E")) +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  xlab("Age") +
  ylab("Force of infection (per susceptible per year)") +
  labs(color    = "Birth cohort era",
       fill     = "Birth cohort era") +
  theme_bw(base_size = 14) +
  theme(legend.position  = "bottom",
        legend.title     = element_text(size = 18),
        legend.text      = element_text(size = 15),
        strip.background = element_rect(fill  = "transparent",
                                        color = "transparent"),
        strip.text       = element_text(size = 12, face = "bold"),
        axis.text.x      = element_text(size = 15),
        axis.text.y      = element_text(size = 15),
      axis.title.y = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         fill  = guide_legend(override.aes = list(alpha = 0.4)))

ggsave("./figures_aug20_2026/foi_curves_by_race_era_study_re.png", #_study_re
       width = 12, height = 8)

## prevalence plot
df_p_curves <- df_p_hat_race_final %>%
  filter(age > 0) %>%  # exclude age 0
  mutate(era = factor(ifelse(birth_cohort_cat == 0, "Pre-1960", "1960+"),
                      levels = c("Pre-1960", "1960+")))

ggplot(df_p_curves,
       aes(x = age, y = Mean, color = era, fill = era, group = era)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.2, color = NA) +
  facet_wrap(~race, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("Pre-1960" = "#336689",
                                "1960+" = "#E65C4E")) +
  scale_fill_manual(values  = c("Pre-1960" = "#336689",
                                "1960+" = "#E65C4E")) +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  xlab("Age") +
  ylab("Force of infection (per susceptible per year)") +
  labs(color    = "Birth cohort era",
       fill     = "Birth cohort era") +
  theme_bw(base_size = 14) +
  theme(legend.position  = "bottom",
        legend.title     = element_text(size = 12),
        legend.text      = element_text(size = 12),
        strip.background = element_rect(fill  = "transparent",
                                        color = "transparent"),
        strip.text       = element_text(size = 12, face = "bold"),
        axis.text.x      = element_text(size = 10),
        axis.text.y      = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         fill  = guide_legend(override.aes = list(alpha = 0.4)))

ggsave("./figures_aug20_2026/p_curves_by_race_era_study_re.png", #_study_re
       width = 12, height = 8)

# *****************************************************************************
#### Computing median age of infection ####
# *****************************************************************************

prev_groups <- list()
foi_groups  <- list()
aoi_groups  <- list()

for (r in race_index) {
  
  alpha_est <- array(NA, dim = c(n_iter, 2))
  gamma_est <- array(NA, dim = c(n_iter, 2))
  age_of_infection_est <- array(NA, dim = c(n_iter, 2))
  p_hat     <- array(NA, dim = c(fin_age + 1, n_iter, 2))
  hp_hat    <- array(NA, dim = c(fin_age + 1, n_iter, 2))
  
  for (iter in 1:n_iter) {
    for (b in 0:1) {
      a <- ilogit(m_alpha0[iter, r] + m_alpha1[iter, r] * b)
      g <- exp(  m_gamma0[iter, r] + m_gamma1[iter, r] * b)
      
      alpha_est[iter, b + 1] <- a
      gamma_est[iter, b + 1] <- g
      age_of_infection_est[iter, b + 1] <- log(2) / g
      
      for (age in 0:fin_age) {
        p_raw <- a * (1 - exp(-g * age))
        p     <- 0.001 + 0.998 * p_raw
        
        p_hat[age + 1, iter, b + 1] <- p_raw
        
        if (age == 0) {
          hp_hat[age + 1, iter, b + 1] <- NA
        } else {
          hp_hat[age + 1, iter, b + 1] <- (a * g * exp(-g * age)) /
                                           (1 - p * (1 - exp(-g * age)) + 1e-10)
        }
      }
    }
  }
  
  dt_foi <- summarize_predictions(hp_hat)
  dt_foi[, `:=`(race = race_labels[r], type = "Force of Infection (\u03bb)")]
  
  # summarize age_of_infection_est directly 
  dt_aoi <- data.frame(
    birth_cohort_cat = 0:1,
    Mean = apply(age_of_infection_est, 2, mean),
    LB   = apply(age_of_infection_est, 2, quantile, 0.025),
    UB   = apply(age_of_infection_est, 2, quantile, 0.975)
  )
  dt_aoi$race <- race_labels[r]
  dt_aoi$type <- "Median age of Infection"
  
  foi_groups[[r]] <- dt_foi
  aoi_groups[[r]] <- dt_aoi
  
  rm(alpha_est, gamma_est, p_hat, hp_hat, dt_foi, dt_aoi)
}

df_aoi_final <- bind_rows(aoi_groups) %>%
  mutate(era = ifelse(birth_cohort_cat == 0, "Pre-1960", "1960+"),
         era = factor(era, levels = c("Pre-1960", "1960+")))

save(df_aoi_final, file='./results/median_age_inf.Rdata')

## plot
ggplot(df_aoi_final,
       aes(x = reorder(race,Mean,decreasing=TRUE), y = Mean, 
           ymin = LB, ymax = UB, group = era, color= race)) +
  #geom_point()+
  geom_pointrange(size = 0.5, position = position_dodge(width = 0.6))+
  geom_errorbar(width = 0.5, size = 1, position = position_dodge(width = 0.6)) +
  #facet_wrap(~era, scales = "fixed")+
  theme_bw(base_size = 14) +
  theme(legend.position  = "bottom", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=14),axis.text.y=element_text(size=15),
        axis.title.x=element_blank(), axis.title.y=element_text(size=15),
      legend.text=element_text(size=10))+ylab("Median age of infection (years)")
  
ggsave("./figures_aug20_2026/median_age_inf_study_re.png") #_study_re

# *****************************************************************************
#### Combined observed points & predicted curves — Prevalence ####
# *****************************************************************************

# observed — age-averaged within era
race_index <- 1:5
race_labels <- c("Hispanic", "NH American Indian", "NH Asian", "NH Black", "NH White")
race_short  <- c("H", "AIAN", "NHA", "NHB", "NHW")
race_titles <- c("Hispanic",
                 "Non-Hispanic American Indian & Alaska Native",
                 "Non-Hispanic Asian",
                 "Non-Hispanic Black",
                 "Non-Hispanic White")

df_p_curves_obs <- data_model %>%
  filter(age > 0) %>%
  mutate(era     = factor(ifelse(birth_cohort_cat == 0, "Pre-1960", "1960+"),
                          levels = c("Pre-1960", "1960+")),
         hp_prev = HP / subgroup_n) %>%
  group_by(race_fac, era, age = round(age)) %>%
  summarise(Mean = mean(hp_prev, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(race=ifelse(race_fac=="Hispanic","Hispanic",
ifelse(race_fac=="Non-Hispanic American Indian","NH American Indian",
ifelse(race_fac=="Non-Hispanic Asian","NH Asian",
ifelse(race_fac=="Non-Hispanic Black", "NH Black",
 "NH White")))))

df_p_curves_obs2 <- data_model_filtered %>%
  filter(age > 0) %>%
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

# predicted
df_p_curves <- df_p_hat_race_final %>%
  filter(age > 0) %>%
  mutate(era = factor(ifelse(birth_cohort_cat == 0, "Pre-1960", "1960+"),
                       levels = c("Pre-1960", "1960+")))
# combined plot
ggplot() +
  geom_ribbon(data = df_p_curves,
              aes(x = age, ymin = LB, ymax = UB,
                  fill = era, group = era),
              alpha = 0.2, color = NA) +
  geom_line(data = df_p_curves,
            aes(x = age, y = Mean, color = era, group = era),
            size = 1) +
  geom_point(data = df_p_curves_obs2,
             aes(x = age, y = Mean, color = era),
             size = 1, alpha = 0.8, shape = 16)+
  geom_errorbar(data = df_p_curves_obs2,
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
ggsave("./figures_aug20_2026/predicted_vs_observed_prev_curves_study_re.png", #_study_re
        width = 12, height = 8)

# *****************************************************************************
#### Visualizations of model convergence ####
# *****************************************************************************

## generate traceplots

generate_traceplot_function <- function(draws,label){
  draws_df <- as.data.frame(draws)
  draws_df$iteration=seq(1,nrow(draws_df))
  draws_df_long <- pivot_longer(draws_df, 
                        cols = -iteration, 
                        names_to = "parameter", 
                        values_to = "value")
  if(label=='alpha intercept') {
   draws_df_long_relabeled=draws_df_long%>%
    mutate(parameter=ifelse(draws_df_long$parameter=="param1","alpha intercept H",
                           ifelse(draws_df_long$parameter=="param2","alpha intercet AIAN",
                           ifelse(draws_df_long$parameter=="param3","alpha intercet NHA",
                           ifelse(draws_df_long$parameter=="param4","alpha intercet NHB","alpha intercept NHW")))))
  }
  else if(label=='alpha slope') {
    draws_df_long_relabeled=draws_df_long%>%
    mutate(parameter=ifelse(draws_df_long$parameter=="param1","alpha slope H",
                           ifelse(draws_df_long$parameter=="param2","alpha slope AIAN",
                           ifelse(draws_df_long$parameter=="param3","alpha slope NHA",
                           ifelse(draws_df_long$parameter=="param4","alpha slope NHB","alpha slope NHW")))))
  }
  else if(label=='gamma intercept') {
     draws_df_long_relabeled=draws_df_long%>%
    mutate(parameter=ifelse(draws_df_long$parameter=="param1","gamma intercept H",
                           ifelse(draws_df_long$parameter=="param2","gamma intercet AIAN",
                           ifelse(draws_df_long$parameter=="param3","gamma intercet NHA",
                           ifelse(draws_df_long$parameter=="param4","gamma intercet NHB","gamma intercept NHW")))))
  }
  else {
    draws_df_long_relabeled=draws_df_long%>%
    mutate(parameter=ifelse(draws_df_long$parameter=="param1","gamma slope H",
                           ifelse(draws_df_long$parameter=="param2","gamma slope AIAN",
                           ifelse(draws_df_long$parameter=="param3","gamma slope NHA",
                           ifelse(draws_df_long$parameter=="param4","gamma slope NHB","gamma slope NHW")))))
  }
  
  ggplot(draws_df_long_relabeled, aes(x = iteration, y = value, color = parameter)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
    theme_minimal() +
    labs(x = "Iteration",
       y = "Estimate") +
  theme(legend.position = "none", 
  axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),strip.text=element_text(size=20))
}

generate_traceplot_function(m_alpha0,"alpha intercept")
ggsave("./figures_aug20_2026/alpha0_traceplot_study_re.png", #_study_re
        width = 12, height = 8)
generate_traceplot_function(m_alpha1,"alpha slope")
ggsave("./figures_aug20_2026/alpha1_traceplot_study_re.png", #_study_re
        width = 12, height = 8)
generate_traceplot_function(m_gamma0,"gamma intercept")
ggsave("./figures_aug20_2026/gamma0_traceplot_study_re.png", #_study_re
        width = 12, height = 8)
generate_traceplot_function(m_gamma1, "gamma slope")
ggsave("./figures_aug20_2026/gamma1_traceplot_study_re.png", #_study_re
        width = 12, height = 8)

## generate density plots
generate_densplot_function <- function(draws,label){
  draws_df <- as.data.frame(draws)
  draws_df$iteration=seq(1,nrow(draws_df))
  draws_df_long <- pivot_longer(draws_df, 
                        cols = -iteration, 
                        names_to = "parameter", 
                        values_to = "value")
  if(label=='alpha intercept') {
   draws_df_long_relabeled=draws_df_long%>%
    mutate(parameter=ifelse(draws_df_long$parameter=="param1","alpha intercept H",
                           ifelse(draws_df_long$parameter=="param2","alpha intercet AIAN",
                           ifelse(draws_df_long$parameter=="param3","alpha intercet NHA",
                           ifelse(draws_df_long$parameter=="param4","alpha intercet NHB","alpha intercept NHW")))))
  }
  else if(label=='alpha slope') {
    draws_df_long_relabeled=draws_df_long%>%
    mutate(parameter=ifelse(draws_df_long$parameter=="param1","alpha slope H",
                           ifelse(draws_df_long$parameter=="param2","alpha slope AIAN",
                           ifelse(draws_df_long$parameter=="param3","alpha slope NHA",
                           ifelse(draws_df_long$parameter=="param4","alpha slope NHB","alpha slope NHW")))))
  }
  else if(label=='gamma intercept') {
     draws_df_long_relabeled=draws_df_long%>%
    mutate(parameter=ifelse(draws_df_long$parameter=="param1","gamma intercept H",
                           ifelse(draws_df_long$parameter=="param2","gamma intercet AIAN",
                           ifelse(draws_df_long$parameter=="param3","gamma intercet NHA",
                           ifelse(draws_df_long$parameter=="param4","gamma intercet NHB","gamma intercept NHW")))))
  }
  else {
    draws_df_long_relabeled=draws_df_long%>%
    mutate(parameter=ifelse(draws_df_long$parameter=="param1","gamma slope H",
                           ifelse(draws_df_long$parameter=="param2","gamma slope AIAN",
                           ifelse(draws_df_long$parameter=="param3","gamma slope NHA",
                           ifelse(draws_df_long$parameter=="param4","gamma slope NHB","gamma slope NHW")))))
  }
  
  ggplot(draws_df_long_relabeled, aes(x = value, fill = parameter)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
    theme_bw() +
    labs(x = "Iteration",
       y = "Estimate") +
  theme(legend.position = "none", 
  axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),strip.text=element_text(size=20))
}

generate_densplot_function(m_alpha0,"alpha intercept")
ggsave("./figures_aug20_2026/alpha0_densplot_study_re.png", #_study_re
        width = 12, height = 8)
generate_densplot_function(m_alpha1,"alpha slope")
ggsave("./figures_aug20_2026/alpha1_densplot_study_re.png", #_study_re
        width = 12, height = 8)
generate_densplot_function(m_gamma0,"gamma intercept")
ggsave("./figures_aug20_2026/gamma0_densplot_study_re.png", #_study_re
        width = 12, height = 8)
generate_densplot_function(m_gamma1, "gamma slope")
ggsave("./figures_aug20_2026/gamma1_densplot_study_re.png", #_study_re
        width = 12, height = 8)

## generate tables for univariate and multivariate convergence statistics

gelman_conv    <- cluster_fit$gelman_conv

gelman_multivariate_PSRF=as.data.frame(gelman_conv$psrf[,1])
colnames(gelman_multivariate_PSRF)="PSRF"
gelman_multivariate_PSRF$param=c("mu_intercept_alpha","mu_slope_alpha",
"mu_intercept_gamma","mu_slope_gamma",
"var_intercept_alpha","var_slope_alpha",
"var_intercept_gamma","var_slope_gamma",
"covar_mu_slope_intercept_alpha","covar_mu_intercept_gamma_alpha",
"covar_mu_intercept_gamma_slope_alpha","covar_slope_gamma_intercept_alpha",
"covar_slope_gamma_alpha","covar_slope_intercept_gamma","overdispersion_parameter","sd_study")

gelman_multivariate_PSRF_table <- gt(cbind.data.frame(parameter=gelman_multivariate_PSRF$param,
PSRF=gelman_multivariate_PSRF$PSRF))
gelman_multivariate_PSRF_table
gtsave(gelman_multivariate_PSRF_table,"./figures_aug20_2026/multivariate_PSRF_study_re.png") #_study_re

gelman_par   <- cluster_fit$gelman_par
gelman_univariate_PSRF=as.data.frame(gelman_par$psrf[,1])
colnames(gelman_univariate_PSRF)="PSRF"
gelman_univariate_PSRF$param=c("log_odds_alpha_H_youngest_age",
"log_odds_alpha_AIAN_youngest_age","log_odds_alpha_NHA_youngest_age",
"log_odds_alpha_NHB_youngest_age","log_odds_alpha_NHW_youngest_age",
"delta_log_odds_alpha_H",
"delta_log_odds_alpha_AIAN","delta_log_odds_alpha_NHA",
"delta_log_odds_alpha_NHB","delta_log_odds_alpha_NHW",
"log_odds_FOI_H_youngest_age",
"log_odds_FOI_AIAN_youngest_age","log_odds_FOI_NHA_youngest_age",
"log_odds_FOI_NHB_youngest_age","log_odds_FOI_NHW_youngest_age",
"delta_log_odds_FOI_H",
"delta_log_odds_FOI_AIAN","delta_log_odds_FOI_NHA",
"delta_log_odds_FOI_NHB","delta_log_odds_FOI_NHW")

gelman_univariate_PSRF_table <- gt(cbind.data.frame(parameter=gelman_univariate_PSRF$param,
PSRF=gelman_univariate_PSRF$PSRF))
gelman_univariate_PSRF_table
gtsave(gelman_univariate_PSRF_table,"./figures_aug20_2026/univariate_PSRF_study_re.png") #_study_re



