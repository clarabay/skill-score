rm(list=ls())

require(arrow)
require(tidyverse)
require(covidHubUtils) # remotes::install_github("reichlab/covidHubUtils")
require(covidData) # remotes::install_github("reichlab/covidData")

### Note that final observed data could change due to data revisions and backfilling 
### compared to the files present in our repository (surveillance data/observed data). 
### To recreate our results, use the surveillance files in the folder of our repo.



### ILI data
# original repository no longer available at the time of publication
# data available here: 'surveillance data/observed data/ilinet_truth.parquet'
#
# library(cdcfluview)
# 
# ilinet_all <- bind_rows(
#   ilinet(region = "national") %>%
#     mutate(location = "US National"),
#   ilinet(region = "hhs") %>%
#     mutate(location = paste("HHS", region))
# )

### influenza hospitalizations
fluh_2223 <- read_csv("https://raw.githubusercontent.com/cdcepi/Flusight-forecast-data/refs/heads/master/data-truth/truth-Incident%20Hospitalizations.csv") %>%
  filter(date <= as.Date('2023-06-10')) %>%
  rename(target_end_date = date)

write_csv(fluh_2223, 'surveillance data/observed data/FluSight-hosp-2022-2023-observed.csv')

fluh_2526 <- read_csv("https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/target-data/target-hospital-admissions.csv") %>%
  filter(as.Date('2023-06-10') < date) %>%
  dplyr::select(-weekly_rate) %>%
  rename(target_end_date = date)

# write_csv(fluh_2526, 'surveillance data/observed data/FluSight-hosp-observed.csv') # older version of data
write_csv(fluh_2526, 'surveillance data/observed data/FluSight-hosp-2025-2026-observed.csv')


### COVID-19 deaths, cases, and daily hospitalizations
df_obs_deaths <- covidHubUtils::load_truth(
  truth_source = "JHU",
  target_variable = "inc death",
  truth_end_date = as.Date('2023-02-01'), 
  temporal_resolution = "weekly",
  locations = hub_locations %>% filter(geo_type == "state") %>% pull(fips)) %>%
  mutate(outcome = 'death') %>%
  dplyr::select(target_end_date, location, location_name, value, outcome)

write_parquet(df_obs_deaths, 'surveillance data/observed data/truth_deaths.parquet')

df_obs_cases <- covidHubUtils::load_truth(
  truth_source = "JHU",
  target_variable = "inc case",
  truth_end_date = as.Date('2022-01-01'), 
  temporal_resolution = "weekly",
  locations = hub_locations %>% pull(fips)) %>%
  mutate(outcome = 'case') %>%
  dplyr::select(target_end_date, location, location_name, value, outcome)

write_parquet(df_obs_cases, 'surveillance data/observed data/truth_cases.parquet')

### COVID-19 forecast hub - DAILY - FINAL as of April 28, 2024
df_obs_hosp_daily <- read_csv('https://media.githubusercontent.com/media/reichlab/covid19-forecast-hub/refs/heads/master/data-truth/truth-Incident%20Hospitalizations.csv') %>%
  mutate(outcome = 'hosp') %>%
  rename(target_end_date = date) 

write_parquet(df_obs_hosp_daily, 'surveillance data/observed data/truth_hosp.parquet')

df_obs_covid <- bind_rows(df_obs_cases, df_obs_deaths, df_obs_hosp_daily) %>%
  mutate(
    obs_value = ifelse(value >= 0, value, NA), 
    log_obs_value = log(obs_value + 1)) %>%
  dplyr::select(location, location_name, outcome, target_end_date, obs_value, log_obs_value)

write_parquet(df_obs_covid, 'surveillance data/observed data/COVID19-observed.parquet')

### COVID-19 weekly hospitalizations
df_obs_hosp_weekly <- read_csv('https://raw.githubusercontent.com/CDCgov/covid19-forecast-hub/refs/heads/main/target-data/covid-hospital-admissions.csv') %>%
  filter(target_end_date <= max(.$target_end_date - 21)) %>%
  mutate(
    obs_value = ifelse(value >= 0, value, NA), 
    log_obs_value = log(obs_value + 1)) %>%
  dplyr::select(location, target_end_date, obs_value, log_obs_value)

write_csv(df_obs_hosp_weekly, 'surveillance data/observed data/COVID19-hosp-observed.csv')



