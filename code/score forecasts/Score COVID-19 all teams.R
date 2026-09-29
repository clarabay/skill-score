rm(list=ls())

require(tidyverse)
require(arrow)

source('code/WIS.R')

locations <- read_csv('dat/locations.csv') %>%
  mutate(location = ifelse(nchar(location) == 1, paste0('0', location), location))

### COVID-19 forecast hub data
df_obs <- read_parquet('surveillance data/observed data/COVID19-observed.parquet')
df_obs_hosp <- read_csv('surveillance data/observed data/COVID19-hosp-observed.csv') %>%
  mutate(outcome = 'hosp') %>%
  left_join(unique(dplyr::select(df_obs, location, location_name)))

df_obs <- bind_rows(df_obs, df_obs_hosp)

# get states
states <- filter(df_obs, outcome == 'hosp', !is.na(location_name)) %>%
  dplyr::select(location, location_name) %>%
  filter(!duplicated(.))#, location != 'US')

### Death, case, and hospitalization forecasts April 2020 to April 2024
all_covid <- open_dataset('dat/covid19-forecast-hub_2020-2024.parquet',
  partitioning = 'team')

### cases
wis_case <- filter(all_covid, type == 'quantile',
  str_detect(target, 'inc'), str_detect(target, 'case')) %>%
  filter(location %in% states$location) %>%
  dplyr::select(team, location, target, quantile, value, forecast_date, target_end_date) %>%
  collect() %>% 
  left_join(filter(df_obs, outcome == 'case') %>% dplyr::select(location, target_end_date, obs_value)) %>%
  group_by(team, location, forecast_date, target_end_date) %>%
  summarize(
    target = first(target),
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log(value + 1), log(obs_value + 1)),
    .groups='drop') %>%
  mutate(
    reference_date = as.Date(ceiling_date(forecast_date, unit="week", week_start=6)),
    reference_date = as.Date(ifelse(forecast_date > (reference_date - 3), reference_date + 7, reference_date)),
    horizon = as.numeric(target_end_date - reference_date)/7,
    #horizon = as.numeric(str_extract(target, '^[0-9]*')),
    #horizon = horizon - 1,
    outcome = 'case'
  ) 
  
write_parquet(wis_case, 'scored forecasts/COVID19-case-ADM1-WIS.parquet')

### daily hosp
wis_hosp <- filter(all_covid, type == 'quantile') %>%
  filter(str_detect(target, 'inc'), str_detect(target, 'hosp')) %>%
  filter(location %in% states$location) %>%
  filter(wday(target_end_date) == 7) %>% # restrict to Saturday outcomes
  dplyr::select(team, location, target, quantile, value, forecast_date, target_end_date) %>%
  collect() %>% 
  left_join(filter(df_obs, outcome == 'hosp') %>% dplyr::select(location, target_end_date, obs_value)) %>%
  group_by(team, location, forecast_date, target_end_date) %>%
  summarize(
    target = first(target),
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log(value + 1), log(obs_value + 1)),
    .groups='drop') %>%
  mutate(
    reference_date = as.Date(ceiling_date(forecast_date, unit="week", week_start=6)),
    # forecasts submitted Wednesday or later are reassigned to the next week
    reference_date = as.Date(ifelse(forecast_date > (reference_date - 3), reference_date + 7, reference_date)),
    horizon = as.numeric(target_end_date - reference_date)/7,
    #horizon = as.numeric(str_extract(target, '^[0-9]*')),
    #horizon = horizon/7 - 1,
    outcome = 'hosp'
  )

write_parquet(wis_hosp, 'scored forecasts/COVID19-hosp-ADM1-WIS.parquet')

### deaths
wis_death <- filter(all_covid, type == 'quantile',
  str_detect(target, 'inc'), str_detect(target, 'death')) %>%
  filter(location %in% states$location) %>%
  dplyr::select(team, location, target, quantile, value, forecast_date, target_end_date) %>%
  collect() %>% 
  left_join(filter(df_obs, outcome == 'death') %>% dplyr::select(location, target_end_date, obs_value)) %>%
  group_by(team, location, forecast_date, target_end_date) %>%
  summarize(
    target = first(target),
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log(value + 1), log(obs_value + 1)),
    .groups='drop') %>%
  mutate(
    reference_date = as.Date(ceiling_date(forecast_date, unit="week", week_start=6)),
    reference_date = as.Date(ifelse(forecast_date > (reference_date - 3), reference_date + 7, reference_date)),
    horizon = as.numeric(target_end_date - reference_date)/7,
    # horizon = as.numeric(str_extract(target, '^[0-9]*')),
    # horizon = horizon - 1,
    outcome = 'death'
  )

write_parquet(wis_death, 'scored forecasts/COVID19-deaths-ADM1-WIS.parquet')


