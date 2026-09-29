rm(list=ls())

require(tidyverse)
require(purrr)
require(arrow)

source('code/forecast_functions.R')

### daily hospitalizations
df_hosp <- read_parquet('surveillance data/observed data/COVID19-observed.parquet') %>%
  filter(outcome == 'hosp') %>%
  filter((as.Date('2020-12-01') - 7) <= target_end_date & target_end_date < (as.Date('2024-06-01') + 7)) %>%
  mutate(obs_value = ifelse(obs_value < 0, NA, obs_value)) %>%
  arrange(target_end_date)
#sum(is.na(df_hosp$obs_value))

quantiles_h <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
  0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)


### version 3 using 3 adjacent data points
hindcast3_h_quant <- tibble()
for (this_loc in unique(df_hosp$location)) {
  for (this_date in unique(df_hosp$target_end_date)) {
    this_date <- as.Date(this_date, origin='1970-01-01')
    if (wday(this_date, label=T) != 'Sat') next()
    this_df <- filter(df_hosp, location == this_loc, 
      target_end_date %in% c(this_date - 1, this_date, this_date + 1))
    if (nrow(this_df) < 3) next()
    if (any(is.na(this_df$obs_value))) next()
    nb_fit <- fit_nbinom(this_df$obs_value, if_all_zeros='add1', 
      weights=c(1, 2, 1), outlier_check=F)
    this_param <- tibble(
      outcome = 'hosp',
      location = this_loc,
      location_name = this_df$location_name[1],
      target_end_date = as.Date(this_date, origin='1970-01-01'),
      nbinom_mu = nb_fit$mu,
      nbinom_size = nb_fit$size,
      outlier_detected = ifelse('outlier_detected' %in% names(nb_fit), nb_fit$outlier_detected, NA)
    )
    this_quant <- this_param %>%
      mutate(qv = map2(nbinom_mu, nbinom_size, ~ nbinom_to_quant(mu=.x, size=.y, quantiles=quantiles_h))) %>%
      unnest(qv) %>%
      dplyr::select(-nbinom_mu, -nbinom_size)
    hindcast3_h_quant <- bind_rows(hindcast3_h_quant, this_quant)
  }
}

write_parquet(hindcast3_h_quant, 'benchmark forecasts/hindcasts/COVID19-hosp-hindcast3-quantiles.parquet')


### cases
df_case <- read_parquet('surveillance data/observed data/COVID19-observed.parquet') %>%
  filter(outcome == 'case') %>%
  filter((as.Date('2020-06-01') - 7) <= target_end_date & target_end_date < (as.Date('2022-01-01') + 7)) %>%
  filter(str_length(location) == 2) %>% # limits to states
  mutate(obs_value = ifelse(obs_value < 0, NA, obs_value)) %>%
  arrange(target_end_date)

quantiles_c <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)

hindcast_c_nb <- tibble()
hindcast_c_quant <- tibble()
for (this_loc in unique(df_case$location)) {
  for (this_week in unique(df_case$target_end_date)) {
    this_df <- filter(df_case, location == this_loc, 
      (this_week - 7) <= target_end_date, target_end_date <= (this_week + 7))
    if (nrow(this_df) < 3) next()
    if (any(is.na(this_df$obs_value))) next()
    nb_fit <- fit_nbinom(this_df$obs_value, if_all_zeros='add1', weights = c(1, 2, 1), outlier_check=F)
    this_param <- tibble(
      outcome = 'case',
      location = this_loc,
      location_name = this_df$location_name[1],
      target_end_date = as.Date(this_week, origin='1970-01-01'),
      nbinom_mu = nb_fit$mu,
      nbinom_size = nb_fit$size,
      outlier_detected = ifelse('outlier_detected' %in% names(nb_fit), nb_fit$outlier_detected, NA)
    )
    this_quant <- this_param %>%
      mutate(qv = map2(nbinom_mu, nbinom_size, ~ nbinom_to_quant(mu=.x, size=.y, quantiles=quantiles_c))) %>%
      unnest(qv) %>%
      dplyr::select(-nbinom_mu, -nbinom_size)
    hindcast_c_nb <- bind_rows(hindcast_c_nb, this_param)
    hindcast_c_quant <- bind_rows(hindcast_c_quant, this_quant)
  }
}

write_parquet(hindcast_c_quant, 'benchmark forecasts/hindcasts/COVID19-case-state-hindcast-quantiles.parquet')  


### deaths
df_death <- read_parquet('surveillance data/observed data/COVID19-observed.parquet') %>%
  filter(outcome == 'death') %>%
  filter((as.Date('2020-04-01') - 7) <= target_end_date & target_end_date < (as.Date('2023-02-01') + 7)) %>%
  mutate(obs_value = ifelse(obs_value < 0, NA, obs_value)) %>%
  arrange(target_end_date) 

quantiles_d <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
  0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

hindcast_d_nb <- tibble()
hindcast_d_quant <- tibble()
for (this_loc in unique(df_death$location)) {
  for (this_week in unique(df_death$target_end_date)) {
    # this_loc = '24'; this_week = as.Date('2021-05-22')
      this_df <- filter(df_death, location == this_loc, 
        (this_week - 7) <= target_end_date, target_end_date <= (this_week + 7))
      if (nrow(this_df) < 3) next()
      if (any(is.na(this_df$obs_value))) next()
      nb_fit <- fit_nbinom(this_df$obs_value, if_all_zeros='add1', weights=c(1, 2, 1), outlier_check=F)
      this_param <- tibble(
        outcome = 'death',
        location = this_loc,
        location_name = this_df$location_name[1],
        target_end_date = as.Date(this_week, origin='1970-01-01'),
        nbinom_mu = nb_fit$mu,
        nbinom_size = nb_fit$size,
        #var_log_obs = var(this_df$log_obs_value),
        outlier_detected = ifelse('outlier_detected' %in% names(nb_fit), nb_fit$outlier_detected, NA)
      )
    this_quant <- this_param %>%
      mutate(qv = map2(nbinom_mu, nbinom_size, ~ nbinom_to_quant(mu=.x, size=.y, quantiles=quantiles_d))) %>%
      unnest(qv) %>%
      dplyr::select(-nbinom_mu, -nbinom_size)
    hindcast_d_nb <- bind_rows(hindcast_d_nb, this_param)
    hindcast_d_quant <- bind_rows(hindcast_d_quant, this_quant)
  }
}

write_parquet(hindcast_d_quant, 'benchmark forecasts/hindcasts/COVID19-death-hindcast-quantiles.parquet')


### weekly hospitalizations 2024-2025
df_hosp_weekly <- read_csv('surveillance data/observed data/COVID19-hosp-observed.csv') %>%
  filter(as.Date('2024-11-09') <= target_end_date) %>%
  mutate(obs_value = ifelse(obs_value < 0, NA, obs_value)) %>%
  arrange(target_end_date)

quantiles_h <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
  0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

hindcast_hweek_quant <- tibble()
for (this_loc in unique(df_hosp_weekly$location)) {
  for (this_date in unique(df_hosp_weekly$target_end_date)) {
    this_df <- filter(df_hosp_weekly, location == this_loc, 
      (this_date - 7) <= target_end_date, target_end_date <= (this_date + 7))
    if (nrow(this_df) < 3) next()
    if (any(is.na(this_df$obs_value))) next()
    nb_fit <- fit_nbinom(this_df$obs_value, if_all_zeros='add1', weights=c(1, 2, 1), outlier_check=F)
    this_param <- tibble(
      outcome = 'hosp',
      location = this_loc,
      target_end_date = as.Date(this_date, origin='1970-01-01'),
      nbinom_mu = nb_fit$mu,
      nbinom_size = nb_fit$size,
      outlier_detected = ifelse('outlier_detected' %in% names(nb_fit), nb_fit$outlier_detected, NA)
    )
    this_quant <- this_param %>%
      mutate(qv = map2(nbinom_mu, nbinom_size, 
        ~ nbinom_to_quant(mu=.x, size=.y, quantiles=quantiles_h))) %>%
      unnest(qv) %>%
      dplyr::select(-nbinom_mu, -nbinom_size)
    hindcast_hweek_quant <- bind_rows(hindcast_hweek_quant, this_quant)
  }
}

write_parquet(hindcast_hweek_quant, 'benchmark forecasts/hindcasts/COVID19-hosp-weekly-hindcast-quantiles.parquet')

