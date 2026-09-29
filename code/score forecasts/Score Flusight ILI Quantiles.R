rm(list=ls())

require(tidyverse)
require(arrow)

source('code/WIS.R')

ilinet_all <- read_parquet('surveillance data/observed data/ilinet_truth.parquet')
ilinet_all <- mutate(ilinet_all, target_end_date = week_start + 6) %>%
  mutate(
    weighted_ili_prop = weighted_ili/100,
    log_weighted_ili_prop = log(weighted_ili_prop + 0.001),
    logit_weighted_ili_prop = boot::logit(weighted_ili_prop),
    ) %>%
  dplyr::select(location, target_end_date, weighted_ili_prop, 
    log_weighted_ili_prop, logit_weighted_ili_prop)

hindcast <- read_parquet('submission format sensitivity analysis/ILI_quantile_forecasts_hindcast_interpolation.parquet') %>%
  filter(str_detect(Model, '23 quantiles')) %>%
  dplyr::select(-Model, -`__index_level_0__`) %>%
  mutate(
    target_end_date = as.Date(target_end_date),
    log_value = log(value + 0.001),
    logit_value = boot::logit(value))

naive <- read_parquet('submission format sensitivity analysis/ILI_quantile_forecasts_naive_interpolation.parquet') %>%
  filter(str_detect(Model, '23 quantiles')) %>%
  dplyr::select(-Model, -`__index_level_0__`) %>%
  mutate(
    forecast_date = as.Date(forecast_date),
    target_end_date = as.Date(target_end_date),
    log_value = log(value + 0.001),
    logit_value = boot::logit(value))

logit_baseline <- read_parquet('submission format sensitivity analysis/ILI_quantile_forecasts_logbaseline_interpolation.parquet') %>%
  filter(str_detect(Model, '23 quantiles')) %>%
  mutate(
    forecast_date = as.Date(forecast_date),
    target_end_date = as.Date(target_end_date),
    log_value = log(value + 0.001),
    logit_value = boot::logit(value)) 

ets_baseline <- read_parquet('submission format sensitivity analysis/ILI_quantile_forecasts_ETS_interpolation.parquet') %>%
  filter(str_detect(Model, '23 quantiles')) %>%
  dplyr::select(-Model, -`__index_level_0__`) %>%
  mutate(
    forecast_date = as.Date(forecast_date),
    target_end_date = as.Date(target_end_date),
    log_value = log(value + 0.001),
    logit_value = boot::logit(value)) 


forecasts <- bind_rows(
    read_csv('dat/ILI team forecasts/ILIforecasts_teams_quantiles_2015-16.csv.gz'),
    read_csv('dat/ILI team forecasts/ILIforecasts_teams_quantiles_2016-17.csv.gz'),
    read_csv('dat/ILI team forecasts/ILIforecasts_teams_quantiles_2017-18.csv.gz'),
    read_csv('dat/ILI team forecasts/ILIforecasts_teams_quantiles_2018-19.csv.gz'),
    read_csv('dat/ILI team forecasts/ILIforecasts_teams_quantiles_2019-20.csv.gz')
  ) %>%
  rename(team = Model) %>%
  filter(team != 'AvgEnsemble') %>%
  bind_rows(
    read_parquet('submission format sensitivity analysis/ILI_quantile_forecasts_Hist-Avg_interpolation.parquet') %>%
      filter(str_detect(Model, '23 quantiles')) %>%
      mutate(
        forecast_date = as.Date(forecast_date),
        target_end_date = as.Date(target_end_date),
        submission_date = as.Date(submission_date),
        team ='Hist-Avg'),
    read_parquet('submission format sensitivity analysis/ILI_quantile_forecasts_interpolation.parquet') %>%
      filter(str_detect(Model, '23 quantiles')) %>%
      mutate(
        forecast_date = as.Date(forecast_date),
        target_end_date = as.Date(target_end_date),
        submission_date = as.Date(submission_date),
        team = 'UnwghtAvg')
  ) %>%
  dplyr::select(-Model, -`__index_level_0__`) %>%
  mutate(logit_value = boot::logit(value)) %>%
  mutate(
    submission_date = as.Date(ifelse(submission_date == '2015-03-17', as.Date('2015-03-16'), submission_date)),
    target_end_date = as.Date(ifelse(submission_date == '2016-12-23', target_end_date + 7, target_end_date)),
    horizon = ifelse(submission_date == '2016-12-23', horizon, horizon))

### fix all forecast dates to ensemble dates
ens_dates <- filter(forecasts, team == 'UnwghtAvg') %>%
  dplyr::select(submission_date, horizon, target_end_date) %>%
  filter(!duplicated(.)) %>%
  rename(forecast_date = submission_date) %>%
  arrange(forecast_date)

forecasts <- dplyr::select(forecasts, -forecast_date) %>%
  left_join(ens_dates) %>%
  arrange(forecast_date)

### WIS
wis_hindcast <- hindcast %>%
  left_join(ilinet_all) %>%
  group_by(location, target_end_date) %>%
  summarize(
    wis_hindcast = weighted_interval_score(quantile, value, weighted_ili_prop),
    logit_wis_hindcast = weighted_interval_score(quantile, logit_value, logit_weighted_ili_prop),
    .groups='drop')

wis_naive <- naive %>%
  left_join(ilinet_all) %>%
  group_by(location, forecast_date, target_end_date) %>%
  summarize(
    horizon = first(horizon),
    wis_naive = weighted_interval_score(quantile, value, weighted_ili_prop),
    logit_wis_naive = weighted_interval_score(quantile, logit_value, logit_weighted_ili_prop),
    .groups='drop')

wis_logit_baseline <- logit_baseline %>%
  left_join(ilinet_all) %>%
  group_by(location, forecast_date, target_end_date) %>%
  summarize(
    wis_logit_baseline = weighted_interval_score(quantile, value, weighted_ili_prop),
    logit_wis_logit_baseline = weighted_interval_score(quantile, logit_value, logit_weighted_ili_prop),
    .groups='drop')

wis_ets_baseline <- ets_baseline %>%
  left_join(ilinet_all) %>%
  group_by(location, forecast_date, target_end_date) %>%
  summarize(
    wis_ets_baseline = weighted_interval_score(quantile, value, weighted_ili_prop),
    logit_wis_ets_baseline = weighted_interval_score(quantile, logit_value, logit_weighted_ili_prop),
    .groups='drop')

wis_forecasts <- forecasts %>%
  left_join(ilinet_all) %>%
  group_by(team, location, forecast_date, target_end_date) %>%
  summarize(
    wis = weighted_interval_score(quantile, value, weighted_ili_prop),
    logit_wis = weighted_interval_score(quantile, logit_value, logit_weighted_ili_prop),
    .groups='drop') %>%
  arrange(forecast_date)

### combined scores
combined <- left_join(wis_forecasts, wis_naive) %>%
  left_join(wis_hindcast) %>%
  left_join(wis_logit_baseline) %>%
  left_join(wis_ets_baseline) %>%
  mutate(  # create and filter to challenge submission dates
    season = ifelse('2015-11-02' <= forecast_date & forecast_date <= '2016-05-16', '2015-2016', NA),
    season = ifelse('2016-11-07' <= forecast_date & forecast_date <= '2017-05-15', '2016-2017', season),
    season = ifelse('2017-11-06' <= forecast_date & forecast_date <= '2018-05-14', '2017-2018', season),
    season = ifelse('2018-10-29' <= forecast_date & forecast_date <= '2019-05-13', '2018-2019', season),
    season = ifelse('2019-10-28' <= forecast_date & forecast_date <= '2020-03-10', '2019-2020', season)
  ) %>%
  filter(!is.na(season)) %>% 
  filter(horizon %in% -1:2) %>% # filter horizon
  group_by(team, location, target_end_date) %>%
  mutate(n_horizon = n()) %>%
  ungroup() %>%
  filter(n_horizon == 4) %>% # filter to only target_end_dates with forecasts at all horizons
  dplyr::select(-n_horizon) %>%
  filter(forecast_date != '2016-12-27') # no ensemble produced

write_parquet(combined, 'scored forecasts/flusight_ILI_combined_WIS_all.parquet')
