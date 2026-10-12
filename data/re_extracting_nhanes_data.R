library(tidyverse)


load(file="data/df_nhanes_88_00_age_group_sex_race.Rdata")
load(file="data/df_nhanes_88_00.Rdata")

df_nhanes_88_00_age_group_sex_race <- df_nhanes_88_00 %>% filter(poll == 1) %>% 
  group_by(age_groups_nh, NHANES, race_f, sex) %>%
  dplyr::summarise(age = mean(age),
                   total = sum(poll),
                   total_w = sum(total_w),
                   total_w_df = sum(total_w_df),
                   birthplaceUS = sum(birthplaceUS),
                   HP = sum(HP),
                   HP_w = sum(HP_w),
                   HP_w_df = sum(HP_w_df),
                   HP_cagA = sum(HP_cagA))  %>% 
  mutate(prev = epitools::binom.exact(x = HP, n = total)[, 3] * 100,
         prev_ci_lb = epitools::binom.exact(x = HP, n = total)[, 4] * 100,
         prev_ci_ub = epitools::binom.exact(x = HP, n = total)[, 5] * 100,
         prev_w = epitools::binom.exact(x = round(HP_w), n = round(total_w))[, 3] * 100,
         prev_ci_lb_w = epitools::binom.exact(x = round(HP_w), n = round(total_w))[, 4] * 100,
         prev_ci_ub_w = epitools::binom.exact(x = round(HP_w), n = round(total_w))[, 5] * 100,
         prev_w_df = epitools::binom.exact(x = round(HP_w_df), n = round(total_w_df))[, 3] * 100,
         prev_ci_lb_w_df = epitools::binom.exact(x = round(HP_w_df), n = round(total_w_df))[, 4] * 100,
         prev_ci_ub_w_df = epitools::binom.exact(x = round(HP_w_df), n = round(total_w_df))[, 5] * 100,
         perc = round((HP/total)*100, digits = 2), #quit
         perc_weights = round((HP_w/total_w)*100, digits = 2),#quit
         perc_cagA = round((HP_cagA/total)*100, digits = 2),#quit
         neg = total - HP, #generate for _w, w_df, cag_w, cag_w_df
         neg_cagA = total - HP_cagA,
         year = case_when(NHANES == "1988-1994" ~ 1991,
                          NHANES == "1999-2000" ~ 2000),
         year_NHANES = year,
         year_period = case_when(year_NHANES == 1991 ~ 0,
                                 year_NHANES == 2000 ~ 9),
         min_age = 0)



df_nhanes_88_00_age_group_sex_race$sex <- factor(df_nhanes_88_00_age_group_sex_race$sex, 
                                                 labels = c("Male",
                                                            "Female"))

df_nhanes_88_00_age_group_sex_race$race_f <- factor(df_nhanes_88_00_age_group_sex_race$race_f, 
                                                    labels = c("Non-Hispanic White",
                                                               "Non-Hispanic Black",
                                                               "Other Hispanic",
                                                               "Mexican-American",
                                                               "Other"))



df_nhanes_88_00_age_sex_race_usa <- df_nhanes_88_00_age_group_sex_race %>% 
                                    mutate(study_id = case_when(NHANES == "1988-1994" ~ "NHANES1994",
                                                                NHANES == "1999-2000" ~ "NHANES2000")) %>% 
   group_by(study_id) %>% mutate(study_n = sum(total),
                    birth_cohort = year - age,
                    race = case_when(race_f == "Other Hispanic"~"Hispanic",
                                    race_f == "Non-Hispanic Black"~"Non-Hispanic Black",
                                     race_f == "Non-Hispanic White"~"Non-Hispanic White",
                                     race_f == "Mexican-American"~"Hispanic",
                                     race_f == "Other"~"Other"),
                    year_pub = case_when(NHANES == "1988-1994" ~ 1995,
                                        NHANES == "1999-2000" ~ 2001),
                    time_period = year,
                    subgroup_n = total) %>% ungroup() %>% 
  dplyr::select(study_id, study_n, subgroup_n, age, time_period, 
         birth_cohort, birthplaceUS, year_pub, year, race_f = race, total, HP) %>% 
group_by(study_id, age, time_period, birth_cohort, year_pub, year, race_f) %>% 
  summarise(study_n = dplyr::first(study_n),
                   across(c(subgroup_n, total, HP), ~ sum(.x)),
                   .groups = "drop") %>% ungroup() %>%
  mutate(prev = HP/total * 100) %>% 
  mutate(prev = epitools::binom.exact(x = round(as.numeric(HP)), n = round(as.numeric(total)))[, 3] * 100,
    prev_ci_lb = epitools::binom.exact(x = round(as.numeric(HP)), n = round(as.numeric(total)))[, 4] * 100,
    prev_ci_ub = epitools::binom.exact(x = round(as.numeric(HP)), n = round(as.numeric(total)))[, 5] * 100)



#
save(df_nhanes_88_00_age_sex_race_usa, file ="data/df_nhanes_88_00_age_sex_race_usa_oct12_2026.Rdata")

load(file = "data/df_hp_usa_w_final.Rdata")

df_hp_usa_w_final <- df_hp_usa_w_final %>%
  filter(!study_id %in% c("NHANES1994", "NHANES2000")) %>%
  bind_rows(df_nhanes_88_00_age_sex_race_usa %>% mutate(version = 1)) %>%
  arrange(year)

# Saved under the same object name (df_hp_usa_w_final) as the file it replaces,
# so every model/figures script's `load(...)` call picks it up as a drop-in
# swap with no other code changes needed.
save(df_hp_usa_w_final, file = "data/df_hp_usa_w_final_unweighted_nhanes.Rdata")
