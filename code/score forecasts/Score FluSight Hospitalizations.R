rm(list=ls())

require(tidyverse)
require(arrow)

source('code/WIS.R')

locations <- read_csv('dat/locations.csv') %>%
  mutate(location = ifelse(nchar(location) == 1, paste0('0', location), location)) %>%
  dplyr::select(-population)

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

hindcast <- read_parquet('benchmark forecasts/hindcasts/FluSight-hosp-hindcast-quantiles.parquet') %>%
  mutate(location_name = ifelse(location_name == "United States", 'US', location_name)) %>%
  mutate(log_value = log(value + 1)) %>%
  arrange(target_end_date)

### team forecasts
# ### 2022-2023
flusight_2223 <- open_dataset('dat/FluSight-hosp-2022-2023-forecast-data.parquet', partitioning='team') %>%
  filter(type != 'point') %>%
  collect() %>%
  mutate(
    log_value = log(value + 1),
    horizon = as.numeric(str_extract(target, '^[1-9]{1}')) - 1,
    reference_date = ceiling_date(forecast_date, unit="week", week_start=6)) %>%
  arrange(target_end_date)
 
# ### 2023 to 2025-26
flusight <- open_dataset('dat/FluSight-hosp-forecast-hub.parquet', partitioning='team') %>%
  filter(output_type == 'quantile', str_detect(target, 'wk inc')) %>%
  collect() %>%
  rename(quantile = output_type_id) %>%
  dplyr::select(-output_type)


flusight <- bind_rows(flusight_2223, flusight) %>%
  mutate(log_value = log(value + 1)) %>%
  arrange(target_end_date)

# ### get WIS
wis_flusight <- flusight %>%
  left_join(df_obs) %>%
  filter(location != '78') %>% # remove USVI, only included to 2022-01-29
  group_by(team, location, reference_date, target_end_date) %>%
  summarize(
    forecast_date = first(forecast_date),
    horizon = first(horizon),
    target = first(target),
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop') %>%
  left_join(locations)

table(wis_flusight$target, str_sub(wis_flusight$target_end_date, 1, 4), is.na(wis_flusight$wis))
# 'MIGHTE-Nsemble' has NA for some quantiles
filter(wis_flusight, str_sub(wis_flusight$target_end_date, 1, 4) == '2023',
  is.na(wis_flusight$wis)) %>%
  dplyr::select(target_end_date, location_name) %>%
  table()
# two weeks for two locations in 2024 missing data
filter(wis_flusight, str_sub(wis_flusight$target_end_date, 1, 4) == '2024',
  is.na(wis_flusight$wis)) %>%
  dplyr::select(target_end_date, location_name) %>%
  table()

wis_flusight <- filter(wis_flusight, !(team == 'MIGHTE-Nsemble' & is.na(wis)))

write_parquet(wis_flusight, 'scored forecasts/FluSight-hosp-WIS.parquet')

# add WIS rolling 'truth'
wis_hindcast <- hindcast %>%
left_join(df_obs) %>%
group_by(location, target_end_date) %>%
summarize(
  wis = weighted_interval_score(quantile, value, obs_value),
  log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
  .groups='drop')
sum(is.na(wis_hindcast$log_wis)) / nrow(wis_hindcast)


### naive forecast

naive <- read_parquet('benchmark forecasts/naive/naive_flu_hosp_h12.parquet') %>%
  mutate(location_name = ifelse(location_name == "United States", 'US', location_name)) %>%
  mutate(log_value = log(value + 1)) %>%
  arrange(target_end_date)

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
    reference_date = as.Date(reference_date)) %>%
  dplyr::select(-target)

sum(is.na(wis_naive$log_wis)) / nrow(wis_naive)


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
 
filter(wis_log_baseline, is.na(wis))


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

### combined scores
wis_forecasts <- read_parquet('scored forecasts/FluSight-hosp-WIS.parquet')

wis_baseline <- filter(wis_forecasts, team == 'FluSight-baseline') %>%
  dplyr::select(location, reference_date, target_end_date, forecast_date, wis, log_wis) %>%
  rename(wis_baseline = wis, log_wis_baseline = log_wis)

combined <- left_join(wis_forecasts, 
  rename(wis_hindcast, wis_hindcast = wis, log_wis_hindcast = log_wis)) %>%
  left_join(rename(wis_naive, wis_naive = wis, log_wis_naive = log_wis)) %>%
  left_join(wis_baseline) %>%
  left_join(dplyr::select(wis_log_baseline, -forecast_date) %>%
      rename(wis_lbaseline = wis, log_wis_lbaseline = log_wis)) %>%
  left_join(dplyr::select(wis_ets_baseline, -forecast_date) %>%
      rename(wis_ets_baseline = wis, log_wis_ets_baseline = log_wis)) %>%

  mutate(  # create and filter to challenge submission dates
    season = ifelse('2022-01-10' <= forecast_date & forecast_date <= '2022-06-20', '2021-2022', NA),
    season = ifelse('2022-10-17' <= forecast_date & forecast_date <= '2023-05-17', '2022-2023', season),
    season = ifelse('2023-10-01' <= reference_date & target_end_date < '2024-05-04', '2023-2024', season),
    season = ifelse('2024-11-20' <= reference_date & reference_date <= '2025-05-31', '2024-2025', season),
    season = ifelse('2025-11-22' <= reference_date & reference_date <= '2026-05-30', '2025-2026', season)
  ) %>%
  filter(!is.na(season)) %>% 
  filter(horizon %in% 0:3) %>% # filter horizon
  group_by(team, location, target_end_date) %>%
  mutate(n_horizon = n()) %>%
  ungroup() %>%
  filter(n_horizon == 4) %>% # filter to only target_end_dates with forecasts at all horizons
  dplyr::select(-n_horizon) %>%
  filter(!(reference_date %in% as.Date(c('2022-01-15', '2022-01-22', '2025-05-31')))) # dates without full data for forecasts


table(filter(combined, str_detect(tolower(team), 'ensemble'))$season)

filter(combined, team == 'FluSight-ensemble', is.na(wis_naive)) %>%
  dplyr::select(season, horizon) %>%
  table()

table(filter(combined, is.na(wis))$target_end_date)

filter(combined, team == 'FluSight-ensemble', is.na(wis_lbaseline)) %>%
    dplyr::select(season, location_name) %>%
  table()

filter(combined, team == 'FluSight-ensemble', is.na(wis_hindcast)) %>%
  dplyr::select(season, target_end_date) %>%
  table()

filter(combined, team == 'FluSight-ensemble', is.na(wis_ets_baseline)) %>%
  dplyr::select(season, target_end_date) %>%
  table()

write_parquet(combined, 'scored forecasts/FluSight-hosp_combined_WIS.parquet')

