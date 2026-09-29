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
  #  filter(location %in% states$location) %>%
  mutate(log_value = log(value + 1)) %>%
  arrange(target_end_date)


naive <- bind_rows(
    read_parquet('benchmark forecasts/naive/naive_covid_cases.parquet'), 
    read_parquet('benchmark forecasts/naive/naive_covid_deaths.parquet'),
    read_parquet('benchmark forecasts/naive/naive_covid_hosp.parquet')
  ) %>%
  mutate(
    log_value = log(value + 1),
    reference_date = ifelse(is.na(reference_date), ceiling_date(forecast_date, unit="week", week_start=6), reference_date),
    reference_date = as.Date(reference_date)) %>%
  arrange(target_end_date)

# ### Death, case, and hospitalization forecasts April 2020 to April 2024
all_covid <- open_dataset('dat/covid19-forecast-hub_2020-2024.parquet',
   partitioning = 'team')

select_teams <- c('COVIDhub-4_week_ensemble', 'COVIDhub-baseline')
select_forecasts <- filter(all_covid, type == 'quantile',
    str_detect(target, 'inc'), target != 'wk inc covid prop ed visits') %>%
    filter(team %in% select_teams) %>%
    filter(location %in% states$location) %>%
    collect() %>%
    mutate(
     log_value = log(value + 1),
     outcome = ifelse(str_detect(target, 'death'), 'death', NA),
     outcome = ifelse(str_detect(target, 'hosp'), 'hosp', outcome),
     outcome = ifelse(str_detect(target, 'case'), 'case', outcome)
    )

table(select_forecasts$team, format(select_forecasts$target_end_date, '%Y'))

### get WIS
wis_subset <- select_forecasts %>%
  left_join(df_obs) %>%
  group_by(team, location, outcome, forecast_date, target_end_date) %>%
  summarize(
    target = first(target),
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop') %>%
    mutate(
    reference_date = as.Date(ceiling_date(forecast_date, unit="week", week_start=6)),
    horizon = as.numeric(target_end_date - reference_date)/7
    # horizon = as.numeric(str_extract(target, '^[0-9]*')),
    # horizon = ifelse(outcome == 'hosp', horizon/7, horizon),
    # horizon = horizon - 1
  )

write_parquet(wis_subset, 'scored forecasts/COVID19-subset-WIS.parquet')


# ### Weekly hospitalization forecasts starting November 2024
covid_hosp_weekly <- open_dataset('dat/covid19-forecast-hub_2024-25.parquet',
  partitioning = 'team')
 
forecasts_hosp <- filter(covid_hosp_weekly, output_type == 'quantile') %>%
  filter(target == 'wk inc covid hosp') %>%
  collect() %>%
  mutate(
    quantile = output_type_id,
    log_value = log(value + 1),
    outcome = 'hosp'
  )

wis_hosp <- forecasts_hosp %>%
  left_join(df_obs) %>%
  group_by(team, location, outcome, reference_date, target_end_date) %>%
  summarize(
    target = first(target),
    horizon = first(horizon),
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop') %>%
  mutate(forecast_date = reference_date - 3)

write_parquet(wis_hosp, 'scored forecasts/COVID19-hosp-WIS.parquet')


# add WIS rolling hindcast
wis_hindcast <- hindcast %>%
  left_join(df_obs) %>%
  group_by(location, outcome, target_end_date) %>%
  summarize(
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop')


### naive forecast
wis_naive <- naive %>%
  left_join(df_obs) %>%
  group_by(location, outcome, reference_date, target_end_date) %>%
  summarize(
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop')

#table(filter(wis_subset, team == 'COVIDhub-4_week_ensemble')$outcome)
table(wis_naive$outcome, format(wis_naive$target_end_date, '%Y'))


### log-baseline forecast
log_baseline <- open_dataset(sources=c(
    'benchmark forecasts/random walk baseline/covid_case_log_baseline_h12.pq', 
    'benchmark forecasts/random walk baseline/covid_dailyhosp_log_baseline_h12.pq', 
    'benchmark forecasts/random walk baseline/covid_death_log_baseline_h12.pq'
  )) %>%
  filter(horizon %in% 1:4 | str_detect(target, 'hosp') & horizon %in% c(5, 12, 19, 26)) %>%
  # includes short daily horizons daily hospitalizations and longer horizons for weekly
  dplyr::select(location, target, forecast_date, target_end_date, quantile, value) %>%
  collect() %>%
  bind_rows(
      read_parquet('benchmark forecasts/random walk baseline/covid_hosp_log_baseline_2025_h12.pq') %>%
        dplyr::select(location, target, reference_date, target_end_date, quantile, value)
    ) %>%
  mutate(
    log_value = log(value + 1),
    #forecast_date = as.Date(forecast_date),
    reference_date = as.Date(reference_date),
    reference_date = as.Date(ifelse(is.na(reference_date), ceiling_date(as.Date(forecast_date), unit="week", week_start=6), reference_date)),
    target_end_date = as.Date(target_end_date),
    #forecast_date = as.Date(ifelse(is.na(forecast_date), reference_date - 3, forecast_date)),
    outcome = ifelse(str_detect(target, 'hosp'), 'hosp', NA),
    outcome = ifelse(str_detect(target, 'death'), 'death', outcome),
    outcome = ifelse(str_detect(target, 'case'), 'case', outcome)
  ) %>%
  arrange(target_end_date)

wis_log_baseline <- log_baseline %>%
  left_join(df_obs) %>%
  group_by(location, forecast_date, reference_date, target_end_date, outcome) %>%
  summarize(
    forecast_date = first(forecast_date),
    target = first(target),
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop') %>%
  dplyr::select(-target, -forecast_date)

### ETS baseline
ets_baseline <- open_dataset(sources=c(
  'benchmark forecasts/trend baseline/covid_case_ets_models.parquet',
  'benchmark forecasts/trend baseline/covid_death_ets_models.parquet',
  'benchmark forecasts/trend baseline/covid_dailyhosp_ets_models.parquet')) %>%
  filter(horizon %in% 1:4 | str_detect(target, 'hosp') & horizon %in% c(5, 12, 19, 26)) %>%
  # includes short daily horizons daily hospitalizations and longer horizons for weekly
  dplyr::select(location, target, forecast_date, target_end_date, quantile, value) %>%
  collect() %>%
    bind_rows(
      read_parquet('benchmark forecasts/trend baseline/covid_weeklyhosp_ets_models.parquet') %>%
        dplyr::select(location, target, reference_date, target_end_date, quantile, value)
    ) %>%
    mutate(
      log_value = log(value + 1),
      #forecast_date = as.Date(forecast_date),
      reference_date = as.Date(reference_date),
      reference_date = as.Date(ifelse(is.na(reference_date), ceiling_date(as.Date(forecast_date), unit="week", week_start=6), reference_date)),
      target_end_date = as.Date(target_end_date),
      #forecast_date = as.Date(ifelse(is.na(forecast_date), reference_date - 3, forecast_date)),
      outcome = ifelse(str_detect(target, 'hosp'), 'hosp', NA),
      outcome = ifelse(str_detect(target, 'death'), 'death', outcome),
      outcome = ifelse(str_detect(target, 'case'), 'case', outcome)
    ) %>%
    arrange(target_end_date)

wis_ets_baseline <- ets_baseline %>%
  left_join(df_obs) %>%
  group_by(location, forecast_date, reference_date, target_end_date, outcome) %>%
  summarize(
    forecast_date = first(forecast_date),
    target = first(target),
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups='drop') %>%
  dplyr::select(-target, -forecast_date)


### combined scores
wis_forecasts <- bind_rows(
    read_parquet('scored forecasts/COVID19-subset-WIS.parquet'), 
    read_parquet('scored forecasts/COVID19-hosp-WIS.parquet')) %>%
  filter(location %in% states$location) %>%
  filter(team %in% c('COVIDhub-4_week_ensemble', 'CovidHub-baseline', 'COVIDhub-baseline',
      'CovidHub-ensemble'))


combined <- left_join(
  wis_forecasts,
    filter(wis_hindcast, location %in% states$location) %>%
      rename(wis_hindcast = wis, log_wis_hindcast = log_wis), 
    by=c('location', 'outcome', 'target_end_date')) %>%
  left_join(
    filter(wis_naive, location %in% states$location) %>%
      rename(wis_naive = wis, log_wis_naive = log_wis),
    by=c('location', 'outcome', 'reference_date', 'target_end_date')) %>%
  left_join(
    filter(wis_forecasts, team %in% c('CovidHub-baseline', 'COVIDhub-baseline'),
        location %in% states$location) %>%
      rename(wis_baseline = wis, log_wis_baseline = log_wis) %>%
      dplyr::select(-team, -target, -horizon, -forecast_date), 
    by=c('location', 'outcome', 'reference_date', 'target_end_date')) %>%
  left_join(
    filter(wis_log_baseline, location %in% states$location) %>%
      rename(wis_lbaseline = wis, log_wis_lbaseline = log_wis), 
    by=c('location', 'outcome', 'reference_date', 'target_end_date')) %>%
  left_join(
    filter(wis_ets_baseline, location %in% states$location) %>%
      rename(wis_ets_baseline = wis, log_wis_ets_baseline = log_wis), 
    by=c('location', 'outcome', 'reference_date', 'target_end_date')) %>%
  mutate(season = format(target_end_date, '%Y')) %>%
  filter(!is.na(wis), !is.na(wis_hindcast)) %>% # this eliminates some) %>%
  left_join(states) %>%
  filter( # filter to challenge submission dates
      (outcome == 'case' & !is.na(forecast_date) & '2020-07-28' <= forecast_date & forecast_date <= '2021-12-21') |
      (outcome == 'death' & !is.na(forecast_date) & '2020-04-27' <= forecast_date & forecast_date <= '2022-12-29') |
      (outcome == 'hosp' & !is.na(forecast_date) & '2021-01-06' <= forecast_date & target_end_date <= '2024-05-04') |
      (outcome == 'hosp' & '2024-11-30' <= reference_date)
  ) %>% 
  filter( # filter to analysis seasons and horizons
    (outcome == 'death' & season %in% c('2020', '2021', '2022')) |
    (outcome == 'case' & season %in% c('2020', '2021')) |
    (outcome == 'hosp' & season %in% c('2020', '2021', '2022', '2023', '2025')),
    horizon %in% 0:3
  ) %>%
  group_by(team, location, outcome, target_end_date) %>%
  mutate(n_horizon = n()) %>%
  ungroup() %>%
  filter(n_horizon == 4) %>% # filter to only target_end_dates with forecasts at all horizons
  dplyr::select(-n_horizon) %>%
  filter(!(location_name %in% c('American Samoa'))) %>% 
  filter(is.na(reference_date) | reference_date != '2025-01-25') %>% # no ensemble produced
  filter(is.na(reference_date) | reference_date < '2025-09-27' | '2025-11-22' <= reference_date) # government shutdown

# 2025: no forecasts for reference data 2025-01-25, though 2 teams made forecasts for this date

write_parquet(combined, 'scored forecasts/COVID19_subset_combined_WIS.parquet')


### all teams
wis_forecasts <- bind_rows(
    read_parquet('scored forecasts/COVID19-case-ADM1-WIS.parquet'),
    read_parquet('scored forecasts/COVID19-hosp-ADM1-WIS.parquet'),
    read_parquet('scored forecasts/COVID19-deaths-ADM1-WIS.parquet'),
    read_parquet('scored forecasts/COVID19-hosp-WIS.parquet')
  )

combined <- left_join(
  wis_forecasts,
    filter(wis_hindcast, location %in% states$location) %>%
      rename(wis_hindcast = wis, log_wis_hindcast = log_wis), 
    by=c('location', 'outcome', 'target_end_date')) %>%
  left_join(
    filter(wis_naive, location %in% states$location) %>%
      rename(wis_naive = wis, log_wis_naive = log_wis),
    by=c('location', 'outcome', 'reference_date', 'target_end_date')) %>%
  left_join(
    filter(wis_forecasts, team %in% c('CovidHub-baseline', 'COVIDhub-baseline'),
        location %in% states$location) %>%
      rename(wis_baseline = wis, log_wis_baseline = log_wis) %>%
      dplyr::select(-team, -target, -horizon, -forecast_date), 
    by=c('location', 'outcome', 'reference_date', 'target_end_date')) %>%
  left_join(
    filter(wis_log_baseline, location %in% states$location) %>%
      rename(wis_lbaseline = wis, log_wis_lbaseline = log_wis), 
    by=c('location', 'outcome', 'reference_date', 'target_end_date')) %>%
  left_join(
    filter(wis_ets_baseline, location %in% states$location) %>%
      rename(wis_ets_baseline = wis, log_wis_ets_baseline = log_wis), 
    by=c('location', 'outcome', 'reference_date', 'target_end_date')) %>%
  mutate(season = format(target_end_date, '%Y')) %>%
  filter(!is.na(wis), !is.na(wis_hindcast)) %>% # this eliminates some) %>%
  left_join(states) %>%
  filter( # filter to challenge submission dates
      (outcome == 'case' & !is.na(forecast_date) & '2020-07-28' <= forecast_date & forecast_date <= '2021-12-21') |
      (outcome == 'death' & !is.na(forecast_date) & '2020-04-27' <= forecast_date & forecast_date <= '2022-12-29') |
      (outcome == 'hosp' & !is.na(forecast_date) & '2021-01-06' <= forecast_date & target_end_date <= '2024-05-04') |
      (outcome == 'hosp' & '2024-11-30' <= reference_date)
  ) %>% 
  filter( # filter to analysis seasons and horizons
    (outcome == 'death' & season %in% c('2020', '2021', '2022')) |
    (outcome == 'case' & season %in% c('2020', '2021')) |
    (outcome == 'hosp' & season %in% c('2020', '2021', '2022', '2023', '2025')),
    horizon %in% 0:3
  ) %>%
  group_by(team, location, outcome, target_end_date) %>%
  mutate(n_horizon = n()) %>%
  ungroup() %>%
  filter(n_horizon == 4) %>% # filter to only target_end_dates with forecasts at all horizons
  dplyr::select(-n_horizon) %>%
  filter(!(location_name %in% c('American Samoa'))) %>% 
  filter(is.na(reference_date) | reference_date != '2025-01-25') %>% # no ensemble produced
  filter(is.na(reference_date) | reference_date < '2025-09-27' | '2025-11-22' <= reference_date) # government shutdown

write_parquet(combined, 'scored forecasts/COVID19_ADM1_combined_WIS.parquet')