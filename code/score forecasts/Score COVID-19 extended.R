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

### load forecasts
hindcast <- bind_rows(
  read_parquet('benchmark forecasts/hindcasts/COVID19-case-state-hindcast-quantiles.parquet'),
  read_parquet('benchmark forecasts/hindcasts/COVID19-hosp-hindcast3-quantiles.parquet'),
  read_parquet('benchmark forecasts/hindcasts/COVID19-hosp-weekly-hindcast-quantiles.parquet'),
  read_parquet('benchmark forecasts/hindcasts/COVID19-death-hindcast-quantiles.parquet')) %>%
  dplyr::select(-location_name) %>%
  filter(location %in% states$location) %>%
  mutate(log_value = log(value + 1)) %>%
  arrange(target_end_date) %>%
  left_join(df_obs) %>%
  group_by(location, outcome, target_end_date) %>%
  summarize(
    #wis = weighted_interval_score(quantile, value, obs_value),
    log_wis_hindcast = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop')

naive <- bind_rows(
    read_parquet('benchmark forecasts/naive/naive_covid_cases_h12.parquet'), 
    read_parquet('benchmark forecasts/naive/naive_covid_deaths_h12.parquet'),
    read_parquet('benchmark forecasts/naive/naive_covid_hosp_h12.parquet')
  ) %>%
  mutate(log_value = log(value + 1)) %>%#,
  arrange(target_end_date) %>%
  left_join(df_obs) %>%
  group_by(location, outcome, horizon, forecast_date, target_end_date) %>%
  summarize(
    reference_date = first(reference_date),
    #wis = weighted_interval_score(quantile, value, obs_value),
    log_wis_naive = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop')

### log-baseline forecast
log_baseline <- open_dataset(sources=c(
    'benchmark forecasts/random walk baseline/covid_case_log_baseline_h12.pq', 
    'benchmark forecasts/random walk baseline/covid_dailyhosp_log_baseline_h12.pq', 
    'benchmark forecasts/random walk baseline/covid_death_log_baseline_h12.pq'
  )) %>%
  # includes short daily horizons daily hospitalizations and longer horizons for weekly
  dplyr::select(location, target, forecast_date, target_end_date, quantile, value) %>%
  collect() %>%
  bind_rows(
      read_parquet('benchmark forecasts/random walk baseline/covid_hosp_log_baseline_2025_h12.pq') %>%
#        rename(reference_date = forecast_date) %>%
        dplyr::select(location, target, reference_date, target_end_date, quantile, value)
    ) %>%
  mutate(
    log_value = log(value + 1),
    forecast_date = as.Date(forecast_date),
    reference_date = as.Date(reference_date),
    target_end_date = as.Date(target_end_date),
    forecast_date = as.Date(ifelse(is.na(forecast_date), reference_date - 3, forecast_date)),
    outcome = ifelse(str_detect(target, 'hosp'), 'hosp', NA),
    outcome = ifelse(str_detect(target, 'death'), 'death', outcome),
    outcome = ifelse(str_detect(target, 'case'), 'case', outcome)
  ) %>%
  arrange(target_end_date) %>%
  left_join(df_obs) %>%
  group_by(location, forecast_date, reference_date, target_end_date, outcome) %>%
  summarize(
    forecast_date = first(forecast_date),
    target = first(target),
    #wis = weighted_interval_score(quantile, value, obs_value),
    log_wis_lbaseline = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop') %>%
  dplyr::select(-target)


### ets forecast
ets_baseline <- open_dataset(sources=c(
  'benchmark forecasts/trend baseline/covid_case_ets_models.parquet',
  'benchmark forecasts/trend baseline/covid_death_ets_models.parquet',
  'benchmark forecasts/trend baseline/covid_dailyhosp_ets_models.parquet')) %>%
  dplyr::select(location, target, forecast_date, target_end_date, quantile, value) %>%
  collect() %>%
  bind_rows(
    read_parquet('benchmark forecasts/trend baseline/covid_weeklyhosp_ets_models.parquet') %>%
      dplyr::select(location, target, reference_date, target_end_date, quantile, value)
  ) %>%
  mutate(
    log_value = log(value + 1),
    forecast_date = as.Date(forecast_date),
    reference_date = as.Date(reference_date),
    target_end_date = as.Date(target_end_date),
    forecast_date = as.Date(ifelse(is.na(forecast_date), reference_date - 3, forecast_date)),
    outcome = ifelse(str_detect(target, 'hosp'), 'hosp', NA),
    outcome = ifelse(str_detect(target, 'death'), 'death', outcome),
    outcome = ifelse(str_detect(target, 'case'), 'case', outcome)
  ) %>%
  arrange(target_end_date) %>%
  left_join(df_obs) %>%
  group_by(location, forecast_date, reference_date, target_end_date, outcome) %>%
  summarize(
    forecast_date = first(forecast_date),
    target = first(target),
    #wis = weighted_interval_score(quantile, value, obs_value),
    log_wis_ets_baseline = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop') %>%
  dplyr::select(-target)

combined <- left_join(naive, log_baseline) %>%
  left_join(hindcast) %>%
  left_join(ets_baseline) %>%
  mutate(season = format(target_end_date, '%Y')) %>%
  left_join(states) %>%
  filter( # filter to challenge submission dates
      (outcome == 'case' & '2020-07-28' <= forecast_date & forecast_date <= '2021-12-21') |
      (outcome == 'death' & '2020-04-27' <= forecast_date & forecast_date <= '2022-12-29') |
      (outcome == 'hosp' & '2021-01-06' <= forecast_date & target_end_date <= '2024-05-04') |
      (outcome == 'hosp' & '2024-11-30' <= reference_date)
  ) %>% 
  filter( # filter to analysis seasons and horizons
    (outcome == 'death' & season %in% c('2020', '2021', '2022')) |
    (outcome == 'case' & season %in% c('2020', '2021')) |
    (outcome == 'hosp' & season %in% c('2020', '2021', '2022', '2023', '2025'))
  ) %>%
  group_by(location, outcome, target_end_date) %>%
  mutate(n_horizon = n()) %>%
  ungroup() %>%
  filter(n_horizon == 12) %>% # filter to only target_end_dates with forecasts at all horizons
  dplyr::select(-n_horizon) %>% 
  filter(is.na(reference_date) | reference_date != '2025-01-25') %>% # no ensemble produced
  filter(is.na(reference_date) | reference_date < '2025-09-27' | '2025-11-22' <= reference_date) %>% # government shutdown
  filter(!(location_name %in% c('American Samoa'))) %>%
  filter(!is.na(log_wis_hindcast))

table(combined$outcome, combined$season)
sum(is.na(combined$log_wis_naive))
sum(is.na(combined$log_wis_hindcast))
sum(is.na(combined$log_wis_lbaseline))
sum(is.na(combined$log_wis_ets_baseline))

filter(combined, is.na(log_wis_lbaseline)) %>%
  dplyr::select(season, outcome) %>%
  table()

write_parquet(combined, 'scored forecasts/COVID19_subset_combined_WIS_h12.parquet')

