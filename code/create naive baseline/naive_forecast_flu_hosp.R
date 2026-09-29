#!/usr/bin/env Rscript
# Generate naive flu hospitalization forecasts for 12 horizons (1-12 wk ahead).
# Same logic as naive_forecast_flu.R but with horizons_to_forecast = 12.
# Output: benchmark forecasts/naive/naive_flu_hosp_h12.parquet
# Run from repo root.
# Use this for skill score calculation vs naive for h1-12.

rm(list = ls())

require(tidyverse)
require(arrow)

source("code/forecast_functions.R")

horizons_to_forecast <- 12

quantiles <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
  0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

name_map <- read_csv("dat/name_map.csv", show_col_types = FALSE)

flu_hosp_versioned <- read_parquet("surveillance data/versioned_data/influenza_hosp.parquet") %>%
  mutate(value = ifelse(value < 0, NA, value)) %>%
  mutate(target_end_date = ceiling_date(time_value, "week", week_start = 7)) %>%
  mutate(target_end_date = target_end_date - 1) %>%
  group_by(as_of, geo_value, target_end_date) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(geo_value) %>%
  filter(target_end_date > "2022-01-01")

last_monday_forecast <- as.Date("2023-05-22")

# ── flu 2022-2023 ───────────────────────────────────────────────────────
df_obs <- read_csv('surveillance data/observed data/FluSight-hosp-2022-2023-observed.csv')  %>%
  filter(as.Date("2021-01-01") <= target_end_date & target_end_date < as.Date("2023-07-01"))

forecast_dates <- unique(df_obs$target_end_date) - 5
forecast_dates <- forecast_dates[as.Date("2022-01-01") <= forecast_dates &
                                   forecast_dates < as.Date("2023-07-01")]

naive_quant <- tibble()
for (this_loc in unique(flu_hosp_versioned$geo_value)) {
  if (which(unique(flu_hosp_versioned$geo_value) == this_loc) %% 10 == 1) {
    message("2022-23: ", this_loc)
  }
  for (this_date in forecast_dates) {
    this_date <- as.Date(this_date, origin = "1970-01-01")
    if (this_date <= last_monday_forecast) {
      this_df <- flu_hosp_versioned %>%
        filter(geo_value == this_loc, as_of == this_date,
               target_end_date <= (this_date - 2), !is.na(value))
    } else {
      this_df <- flu_hosp_versioned %>%
        filter(geo_value == this_loc, as_of == this_date,
               target_end_date <= (this_date - 4), !is.na(value))
    }
    if (nrow(this_df) < 3) next()
    nb_fit <- fit_nbinom(this_df$value, if_all_zeros = "add1")
    this_quant <- tibble()
    for (h in 1:horizons_to_forecast) {
      row_param <- tibble(
        location = name_map$location_number[which(name_map$geo_value == this_loc)],
        location_name = name_map$location_name[which(name_map$geo_value == this_loc)],
        forecast_date = this_date,
        reference_date = as.Date(NA),
        target = paste0(h, " wk ahead inc flu hosp"),
        target_end_date = as.Date(this_date, origin = "1970-01-01") + 5 + (h - 1) * 7,
        nbinom_mu = nb_fit$mu,
        nbinom_size = nb_fit$size
      )
      row_quant <- row_param %>%
        mutate(qv = map2(nbinom_mu, nbinom_size,
                        ~ nbinom_to_quant(mu = .x, size = .y, quantiles = quantiles))) %>%
        unnest(qv) %>%
        dplyr::select(-nbinom_mu, -nbinom_size)
      this_quant <- bind_rows(this_quant, row_quant)
    }
    naive_quant <- bind_rows(naive_quant, this_quant)
  }
}

# ── flu 2023-2024 ───────────────────────────────────────────────────────
df_hosp_2324 <- read_csv('surveillance data/observed data/FluSight-hosp-2025-2026-observed.csv') %>%
  mutate(
    target_end_date = as.Date(target_end_date),
    obs_value = ifelse(value < 0, NA, value)
  ) %>%
  arrange(target_end_date)

forecast_dates_2324 <- unique(df_hosp_2324$target_end_date) - 3
forecast_dates_2324 <- forecast_dates_2324[as.Date("2023-10-01") < forecast_dates_2324 &
                                             forecast_dates_2324 < as.Date("2024-06-01")]

naive_2324_quant <- tibble()
for (this_loc in unique(flu_hosp_versioned$geo_value)) {
  if (which(unique(flu_hosp_versioned$geo_value) == this_loc) %% 10 == 1) {
    message("2023-24: ", this_loc)
  }
  for (this_date in forecast_dates_2324) {
    this_date <- as.Date(this_date, origin = "1970-01-01")
    this_df <- flu_hosp_versioned %>%
      filter(geo_value == this_loc, as_of == this_date,
             target_end_date < this_date, !is.na(value))
    if (nrow(this_df) < 3) next()
    nb_fit <- fit_nbinom(this_df$value, if_all_zeros = "add1")
    for (h in 1:horizons_to_forecast) {
      row_param <- tibble(
        location = name_map$location_number[which(name_map$geo_value == this_loc)],
        location_name = name_map$location_name[which(name_map$geo_value == this_loc)],
        forecast_date = this_date,
        reference_date = as.Date(this_date, origin = "1970-01-01") + 3,
        target = paste0(h, " wk ahead inc flu hosp"),
        target_end_date = as.Date(this_date, origin = "1970-01-01") + 3 + (h - 1) * 7,
        nbinom_mu = nb_fit$mu,
        nbinom_size = nb_fit$size
      )
      row_quant <- row_param %>%
        mutate(qv = map2(nbinom_mu, nbinom_size,
                        ~ nbinom_to_quant(mu = .x, size = .y, quantiles = quantiles))) %>%
        unnest(qv) %>%
        dplyr::select(-nbinom_mu, -nbinom_size)
      naive_2324_quant <- bind_rows(naive_2324_quant, row_quant)
    }
  }
}

# ── flu 2024-2025 ───────────────────────────────────────────────────────
df_hosp_2425 <- read_parquet("surveillance data/versioned_data/influenza_hosp_24-26.parquet") %>%
  mutate(
    target_end_date = as.Date(time_value) + 6,
    obs_value = ifelse(value < 0, NA, value)
  ) %>%
  filter(target_end_date <= "2024-05-04" | target_end_date >= "2024-11-09") %>%
  arrange(target_end_date)

forecast_dates_2425 <- unique(df_hosp_2425$as_of)
forecast_dates_2425 <- forecast_dates_2425[as.Date("2024-11-20") <= forecast_dates_2425]

naive_2425_quant <- tibble()
for (this_loc in unique(df_hosp_2425$geo_value)) {
  if (which(unique(df_hosp_2425$geo_value) == this_loc) %% 10 == 1) {
    message("2024-25: ", this_loc)
  }
  for (this_date in forecast_dates_2425) {
    this_date <- as.Date(this_date, origin = "1970-01-01")
    this_df <- df_hosp_2425 %>%
      filter(geo_value == this_loc, as_of == this_date,
             target_end_date < this_date, !is.na(value))
    if (nrow(this_df) < 3) next()
    nb_fit <- fit_nbinom(this_df$value, if_all_zeros = "add1")
    for (h in 1:horizons_to_forecast) {
      row_param <- tibble(
        location = name_map$location_number[which(name_map$geo_value == this_loc)],
        location_name = name_map$location_name[which(name_map$geo_value == this_loc)],
        forecast_date = this_date,
        reference_date = as.Date(this_date, origin = "1970-01-01") + 3,
        target = paste0(h, " wk ahead inc flu hosp"),
        target_end_date = as.Date(this_date, origin = "1970-01-01") + 3 + (h - 1) * 7,
        nbinom_mu = nb_fit$mu,
        nbinom_size = nb_fit$size
      )
      row_quant <- row_param %>%
        mutate(qv = map2(nbinom_mu, nbinom_size,
                        ~ nbinom_to_quant(mu = .x, size = .y, quantiles = quantiles))) %>%
        unnest(qv) %>%
        dplyr::select(-nbinom_mu, -nbinom_size)
      naive_2425_quant <- bind_rows(naive_2425_quant, row_quant)
    }
  }
}

# ── Combine and save ─────────────────────────────────────────────────────
naive_flu_hosp_h12 <- bind_rows(naive_quant, naive_2324_quant, naive_2425_quant)

write_parquet(naive_flu_hosp_h12, "benchmark forecasts/naive/naive_flu_hosp_h12.parquet")
message("Saved benchmark forecasts/naive/naive_flu_hosp_h12.parquet")
message("Horizons: 1-", horizons_to_forecast, " wk ahead; rows: ", nrow(naive_flu_hosp_h12))
