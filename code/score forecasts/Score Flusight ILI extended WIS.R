rm(list=ls())

require(tidyverse)
require(arrow)

source('code/WIS.R')

ilinet_all <- read_parquet('surveillance data/observed data/ilinet_truth.parquet')
ilinet_all <- mutate(ilinet_all, target_end_date = week_start + 6) %>%
  mutate(
    weighted_ili_prop = weighted_ili/100,
    logit_weighted_ili_prop = boot::logit(weighted_ili_prop),
  ) %>%
  dplyr::select(location, target_end_date, weighted_ili_prop, 
    logit_weighted_ili_prop)

hindcast <- read_parquet('submission format sensitivity analysis/ILI_quantile_forecasts_hindcast_interpolation.parquet') %>%
  filter(str_detect(Model, '23 quantiles')) %>%
  dplyr::select(-Model, -`__index_level_0__`) %>%
  mutate(
    target_end_date = as.Date(target_end_date),
    logit_value = boot::logit(value)) %>%
  left_join(ilinet_all) %>%
  group_by(location, target_end_date) %>%
  summarize(
    wis_hindcast = weighted_interval_score(quantile, value, weighted_ili_prop),
    logit_wis_hindcast = weighted_interval_score(quantile, logit_value, logit_weighted_ili_prop),
    .groups='drop')

naive <- read_parquet('submission format sensitivity analysis/ILI_quantile_forecasts_naive_interpolation.parquet') %>%
  filter(str_detect(Model, '23 quantiles')) %>%
  dplyr::select(-Model, -`__index_level_0__`) %>%
  mutate(
    forecast_date = as.Date(forecast_date),
    target_end_date = as.Date(target_end_date),
    logit_value = boot::logit(value)) %>%
  left_join(ilinet_all) %>%
  group_by(location, forecast_date, target_end_date) %>%
  summarize(
    horizon = first(horizon),
    wis_naive = weighted_interval_score(quantile, value, weighted_ili_prop),
    logit_wis_naive = weighted_interval_score(quantile, logit_value, logit_weighted_ili_prop),
    .groups='drop')

logit_baseline <- read_parquet('submission format sensitivity analysis/ILI_quantile_forecasts_logbaseline_interpolation.parquet') %>%
  filter(str_detect(Model, '23 quantiles')) %>%
  mutate(
    forecast_date = as.Date(forecast_date),
    target_end_date = as.Date(target_end_date),
    logit_value = boot::logit(value)) %>%
  left_join(ilinet_all) %>%
  group_by(location, forecast_date, target_end_date) %>%
  summarize(
    logit_wis_logit_baseline = weighted_interval_score(quantile, logit_value, logit_weighted_ili_prop),
    .groups='drop')

ets_baseline <- read_parquet('submission format sensitivity analysis/ILI_quantile_forecasts_ETS_interpolation.parquet') %>%
  filter(str_detect(Model, '23 quantiles')) %>%
  dplyr::select(-Model, -`__index_level_0__`) %>%
  mutate(
    forecast_date = as.Date(forecast_date),
    target_end_date = as.Date(target_end_date),
    logit_value = boot::logit(value)) %>%
  left_join(ilinet_all) %>%
  group_by(location, forecast_date, target_end_date) %>%
  summarize(
    logit_wis_ets_baseline = weighted_interval_score(quantile, logit_value, logit_weighted_ili_prop),
    .groups='drop')

seasonal_baseline <- read_parquet('submission format sensitivity analysis/ILI_quantile_forecasts_Hist-Avg_interpolation.parquet') %>%
    filter(str_detect(Model, '23 quantiles')) %>%
    mutate(
      forecast_date = as.Date(forecast_date),
      target_end_date = as.Date(target_end_date),
      logit_value = boot::logit(value)) %>%
    left_join(ilinet_all) %>%
    group_by(location, forecast_date, target_end_date) %>%
    summarize(
      wis_seas_baseline = weighted_interval_score(quantile, value, weighted_ili_prop),
      logit_wis_seas_baseline = weighted_interval_score(quantile, logit_value, logit_weighted_ili_prop),
      .groups='drop') %>%
  dplyr::select(-forecast_date) %>%
  filter(!duplicated(dplyr::select(., location, target_end_date)))

### combined scores
combined <- left_join(logit_baseline, naive) %>%
  left_join(ets_baseline) %>%
  left_join(hindcast) %>%
  left_join(seasonal_baseline) %>%
  mutate(  # create and filter to challenge submission dates
    season = ifelse('2014-10-20' <= forecast_date & forecast_date <= '2015-05-25', '2014-2015', NA),
    season = ifelse('2015-11-02' <= forecast_date & forecast_date <= '2016-05-16', '2015-2016', season),
    season = ifelse('2016-11-07' <= forecast_date & forecast_date <= '2017-05-15', '2016-2017', season),
    season = ifelse('2017-11-06' <= forecast_date & forecast_date <= '2018-05-14', '2017-2018', season),
    season = ifelse('2018-10-29' <= forecast_date & forecast_date <= '2019-05-13', '2018-2019', season),
    season = ifelse('2019-10-28' <= forecast_date & forecast_date <= '2020-03-10', '2019-2020', season)
  ) %>%
  filter(!is.na(season)) %>% 
#  filter(horizon %in% 0:3) %>% # filter horizon
  group_by(location, target_end_date) %>%
  mutate(n_horizon = n()) %>%
  ungroup() %>%
  filter(n_horizon == 12) %>% # filter to only target_end_dates with forecasts at all horizons
  dplyr::select(-n_horizon)

table(combined$season)
sum(is.na(combined$logit_wis_naive))
sum(is.na(combined$logit_wis_logit_baseline))
sum(is.na(combined$logit_wis_ets_baseline))
sum(is.na(combined$logit_wis_hindcast))

write_csv(combined, 'scored forecasts/flusight_combined_WIS_h12.csv')

