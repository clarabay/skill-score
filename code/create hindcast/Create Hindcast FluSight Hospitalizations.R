rm(list=ls())

require(tidyverse)
require(fitdistrplus)
require(purrr)
require(arrow)

source('code/forecast_functions.R')

df_obs <- bind_rows(
  filter(read_csv('surveillance data/observed data/FluSight-hosp-2022-2023-observed.csv'), target_end_date <= '2023-06-10'),
  filter(read_csv('surveillance data/observed data/FluSight-hosp-2025-2026-observed.csv'), '2023-06-10' < target_end_date)
) %>%
  mutate(
    location_name = ifelse(location_name == 'US', 'United States', location_name),
    obs_value = ifelse(value < 0, NA, value)) %>% 
  dplyr::select(-value) %>%
  filter(target_end_date <= '2024-05-04' | '2024-11-09' <= target_end_date) %>%
  arrange(target_end_date)

fcast_quantiles <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 
  0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

hindcast_param <- tibble()
hindcast_quant <- tibble()
for (this_loc in unique(df_obs$location)) {
  for (this_date in unique(df_obs$target_end_date)) {
    this_date <- as.Date(this_date, origin='1970-01-01')
    this_df <- filter(df_obs, location == this_loc, 
      (this_date - 7) <= target_end_date, target_end_date <= (this_date + 7))
    if (nrow(this_df) < 3) next()
    if (any(is.na(this_df$obs_value))) next()
    nb_fit <- fit_nbinom(this_df$obs_value, if_all_zeros='add1', weights=c(1, 2, 1), outlier_check=F)
    this_param <- tibble(
      location = this_loc,
      location_name = this_df$location_name[1],
      target_end_date = as.Date(this_date, origin='1970-01-01'),
      nbinom_mu = nb_fit$mu,
      nbinom_size = nb_fit$size,
      outlier_detected = ifelse('outlier_detected' %in% names(nb_fit), nb_fit$outlier_detected, NA)
    )
    this_quant <- this_param %>%
      mutate(qv = map2(nbinom_mu, nbinom_size, 
        ~ nbinom_to_quant(mu=.x, size=.y, quantiles=fcast_quantiles))) %>%
      unnest(qv) %>%
      dplyr::select(-nbinom_mu, -nbinom_size)
    if (any(is.na(this_quant$value))) stop()
    hindcast_param <- bind_rows(hindcast_param, this_param)
    hindcast_quant <- bind_rows(hindcast_quant, this_quant)
  }
}

write_parquet(hindcast_quant, 'benchmark forecasts/hindcasts/FluSight-hosp-hindcast-quantiles.parquet')

