rm(list=ls())

require(tidyverse)
require(arrow)

source('code/WIS.R')

locations <- read_csv('dat/locations.csv') %>%
  mutate(location = ifelse(nchar(location) == 1, paste0('0', location), location))

df_obs <- bind_rows(
  filter(read_csv('surveillance data/observed data/FluSight-hosp-2022-2023-observed.csv'), target_end_date <= '2023-06-10'),
  filter(read_csv('surveillance data/observed data/FluSight-hosp-2025-2026-observed.csv'), '2023-06-10' < target_end_date)
) %>%
  mutate(
    #location_name = ifelse(location_name == 'US', 'United States', location_name),
    obs_value = ifelse(value < 0, NA, value),
    log_obs_value = log(obs_value + 1)) %>% 
  dplyr::select(-value, -location_name) %>%
  arrange(target_end_date)

### hindcast
hindcast <- read_parquet('benchmark forecasts/hindcasts/FluSight-hosp-hindcast-quantiles.parquet') %>%
  mutate(location_name = ifelse(location_name == "United States", 'US', location_name)) %>%
  mutate(log_value = log(value + 1)) %>%
  arrange(target_end_date)

wis_hindcast <- hindcast %>%
  left_join(df_obs) %>%
  group_by(location, target_end_date) %>%
  summarize(
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop')

### naive
naive <- read_parquet('benchmark forecasts/naive/naive_flu_hosp_h12.parquet') %>%
  mutate(location_name = ifelse(location_name == "United States", 'US', location_name)) %>%
  mutate(log_value = log(value + 1)) %>%
  arrange(target_end_date) %>%
  dplyr::select(-location_name)

wis_naive <- naive %>%
  left_join(df_obs) %>%
  group_by(location, forecast_date, reference_date, target_end_date) %>%
  summarize(
    target = first(target),
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop') %>%
  mutate(
    reference_date = ifelse(is.na(reference_date), ceiling_date(forecast_date, unit="week", week_start=6), reference_date),
    reference_date = as.Date(reference_date),
    horizon = (target_end_date - reference_date)/7) %>%
  dplyr::select(-target)

### log-baseline forecast
log_baseline <- read_parquet('benchmark forecasts/random walk baseline/flu_hosp_log_baseline_h12.pq') %>%
  mutate(
    log_value = log(value + 1),
    forecast_date = as.Date(forecast_date),
    reference_date = as.Date(reference_date),
    reference_date = ifelse(is.na(reference_date), ceiling_date(forecast_date, unit="week", week_start=6), reference_date),
    reference_date = as.Date(reference_date)) %>%
  arrange(target_end_date)

wis_log_baseline <- log_baseline %>%
  left_join(df_obs) %>%
  group_by(location, reference_date, target_end_date) %>%
  summarize(
    forecast_date = first(forecast_date),
    target = first(target),
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop') %>%
  dplyr::select(-target)


### ets baseline
ets_baseline <- read_parquet('benchmark forecasts/trend baseline/flu_weeklyhosp_ets_models.parquet') %>%
  mutate(
    log_value = log(value + 1),
    forecast_date = as.Date(forecast_date),
    reference_date = ceiling_date(forecast_date, unit="week", week_start=6)
  ) %>%
  arrange(target_end_date) %>%
  dplyr::select(-location_name)

wis_ets_baseline <- ets_baseline %>%
  left_join(df_obs) %>%
  group_by(location, reference_date, target_end_date) %>%
  summarize(
    forecast_date = first(forecast_date),
    target = first(target),
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop') %>%
  dplyr::select(-target)

### combine
combined <- left_join(
    rename(wis_naive, wis_naive = wis, log_wis_naive = log_wis), 
    rename(wis_hindcast, wis_hindcast = wis, log_wis_hindcast = log_wis)) %>%
  left_join(dplyr::select(wis_log_baseline, -forecast_date) %>%
      rename(wis_lbaseline = wis, log_wis_lbaseline = log_wis)) %>%
  left_join(dplyr::select(wis_ets_baseline, -forecast_date) %>%
      rename(wis_ets_baseline = wis, log_wis_ets_baseline = log_wis)) %>%
  left_join(locations) %>%
  mutate(  # create and filter to challenge submission dates
    season = ifelse('2022-01-10' <= forecast_date & forecast_date <= '2022-06-20', '2021-2022', NA),
    season = ifelse('2022-10-17' <= forecast_date & forecast_date <= '2023-05-17', '2022-2023', season),
    season = ifelse('2023-10-01' <= reference_date & target_end_date < '2024-05-04', '2023-2024', season),
    season = ifelse('2024-11-20' <= reference_date & reference_date <= '2025-05-31', '2024-2025', season),
    season = ifelse('2025-11-22' <= reference_date & reference_date <= '2026-05-30', '2025-2026', season)
  ) %>%
  filter(!is.na(season)) %>% 
  group_by(location, target_end_date) %>%
  mutate(n_horizon = n()) %>%
  ungroup() %>%
  filter(n_horizon == 12) %>% # filter to only target_end_dates with forecasts at all horizons
  dplyr::select(-n_horizon) %>%
  filter(!(reference_date %in% as.Date(c('2022-01-15', '2022-01-22', '2025-05-31')))) %>% # dates without full data for forecasts
  filter(!(reference_date %in% c('2023-10-07', '2025-01-25'))) %>% # no ensemble and thus no log_baseline generated
  filter(!(location_name %in% c('American Samoa', 'Guam', 'Northern Mariana Islands', 
    'Virgin Islands', 'Puerto Rico')))
  
table(combined$season)

sum(is.na(combined$log_wis_naive))
sum(is.na(combined$log_wis_hindcast))
sum(is.na(combined$log_wis_lbaseline))
sum(is.na(combined$log_wis_ets_baseline))

# filter(combined, is.na(wis_lbaseline)) %>%
#   dplyr::select(location_name, season) %>%
#   table()

write_parquet(combined, 'scored forecasts/FluSight-hosp_combined_WIS_h12.parquet')




