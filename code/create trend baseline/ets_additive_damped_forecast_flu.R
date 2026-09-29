#!/usr/bin/env Rscript
# Additive damped ETS (AAN, damped=TRUE) for FluSight hospitalization forecasts.
# Four outputs: (STL × log) × {on, off}.
# Run from repository root:  Rscript code/ets_additive_damped_forecast_flu.R

rm(list = ls())

args      <- commandArgs(trailingOnly = TRUE)
only_log  <- "--only-log" %in% args

require(tidyverse)
require(arrow)
require(forecast)

source("code/create trend baseline/ets_forecast_flu_helpers.R")

name_map <- read.csv("dat/name_map.csv")

df_obs <- read_csv('surveillance data/observed data/FluSight-hosp-2022-2023-observed.csv')  %>%
  filter(as.Date("2021-01-01", origin='1970-01-01') <= target_end_date & target_end_date < as.Date("2023-07-01", origin='1970-01-01')) %>%
  mutate(obs_value = ifelse(value < 0, NA, value)) %>%
  arrange(target_end_date)

forecast_dates <- unique(df_obs$target_end_date) - 5
forecast_dates <- forecast_dates[
  as.Date("2022-10-01", origin='1970-01-01') <= forecast_dates & forecast_dates < as.Date("2023-07-01", origin='1970-01-01')
]

flu_hosp_versioned <- read_parquet("surveillance data/versioned_data/influenza_hosp.parquet") %>%
  mutate(value = ifelse(value < 0, NA, value)) %>%
  mutate(target_end_date = ceiling_date(time_value, "week", week_start = 7)) %>%
  mutate(target_end_date = target_end_date - 1) %>%
  group_by(as_of, geo_value, target_end_date) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(geo_value) %>%
  filter(target_end_date > "2022-01-01")

flu_hosp_versioned_2122 <- read_parquet("surveillance data/versioned_data/influenza_hosp.parquet") %>%
  mutate(value = ifelse(value < 0, NA, value)) %>%
  mutate(target_end_date = ceiling_date(time_value, "week", week_start = 7)) %>%
  mutate(target_end_date = target_end_date - 1) %>%
  group_by(as_of, geo_value, target_end_date) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(geo_value) %>%
  filter(target_end_date > "2021-08-01")

last_monday_forecast <- as.Date("2023-05-22")
forecast_dates_2122 <- seq(as.Date("2022-01-10"), as.Date("2022-06-20"), by = "week")

df_hosp_2324 <- read_csv('surveillance data/observed data/FluSight-hosp-2025-2026-observed.csv') %>%
  mutate(
    target_end_date = as.Date(target_end_date),
    obs_value = ifelse(value < 0, NA, value)
  ) %>%
  arrange(target_end_date)

forecast_dates_2324 <- unique(df_hosp_2324$target_end_date) - 3
forecast_dates_2324 <- forecast_dates_2324[
  as.Date("2023-10-01") < forecast_dates_2324 & forecast_dates_2324 < as.Date("2024-06-01")
]

df_hosp_2425 <- read_parquet("surveillance data/versioned_data/influenza_hosp_24-26.parquet") %>%
  mutate(target_end_date = as.Date(time_value) + 6) %>%
  arrange(target_end_date)

forecast_dates_2425 <- unique(df_hosp_2425$as_of)
forecast_dates_2425 <- forecast_dates_2425[
  as.Date("2024-11-20") <= forecast_dates_2425
]

# --- core: one location-date forecast ---
forecast_additive_one <- function(y_raw, meta, use_stl, use_log) {
  if (length(y_raw) < 4L) return(NULL)
  stl_ok <- isTRUE(use_stl) && length(y_raw) >= 2L * STL_PERIOD

  if (use_log) {
    stl_log <- if (stl_ok) stl_prefilter_log1p(y_raw, STL_PERIOD, horizons_to_forecast) else NULL
    x <- log(y_raw + 1)
    if (!is.null(stl_log)) x <- stl_log$x_deseasoned
    adj <- if (!is.null(stl_log)) stl_log$S_future else rep(0, horizons_to_forecast)
    if (all(x == 0)) x <- c(x, log(1 + 1))
    if (any(!is.finite(x))) return(NULL)
    m <- tryCatch(ets(x, model = "AAN", damped = TRUE), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    pred <- tryCatch(
      forecast(m, h = horizons_to_forecast, level = forecast_levels),
      error = function(e) NULL
    )
    if (is.null(pred)) return(NULL)
    pred_to_quant(pred, meta, log_scale = TRUE, forecast_adjustment = adj)
  } else {
    stl_lev <- if (stl_ok) stl_prefilter_level_add(y_raw, STL_PERIOD, horizons_to_forecast) else NULL
    x <- y_raw + 1
    if (!is.null(stl_lev)) x <- stl_lev$x_deseasoned
    adj <- if (!is.null(stl_lev)) stl_lev$S_future else rep(0, horizons_to_forecast)
    if (any(!is.finite(x))) return(NULL)
    m <- tryCatch(ets(x, model = "AAN", damped = TRUE), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    pred <- tryCatch(
      forecast(m, h = horizons_to_forecast, level = forecast_levels),
      error = function(e) NULL
    )
    if (is.null(pred)) return(NULL)
    pred_to_quant(pred, meta, log_scale = FALSE, forecast_adjustment = adj)
  }
}

run_season_2122 <- function(use_stl, use_log) {
  out <- tibble()
  for (this_loc in unique(flu_hosp_versioned_2122$geo_value)) {
    print(this_loc)
    for (this_date in forecast_dates_2122) {
      this_date <- as.Date(this_date, origin='1970-01-01')
      this_df <- flu_hosp_versioned_2122 %>%
        filter(
          geo_value == this_loc, as_of == this_date,
          target_end_date <= (this_date - 2), !is.na(value)
        )
      if (nrow(this_df) < 4) next()
      meta <- list(
        location = name_map$location_number[which(name_map$geo_value == this_loc)],
        location_name = name_map$location_name[which(name_map$geo_value == this_loc)],
        forecast_date = as.Date(this_date, origin = "1970-01-01"),
        reference_date = as.Date(NA),
        target_end_date_base = as.Date(this_date, origin = "1970-01-01") + 5
      )
      q <- forecast_additive_one(this_df$value, meta, use_stl, use_log)
      if (!is.null(q)) out <- bind_rows(out, q)
    }
  }
  out
}

run_season_2223 <- function(use_stl, use_log) {
  out <- tibble()
  for (this_loc in unique(flu_hosp_versioned$geo_value)) {
    print(this_loc)
    for (this_date in forecast_dates) {
      this_date <- as.Date(this_date, origin='1970-01-01')
      if (this_date <= last_monday_forecast) {
        this_df <- flu_hosp_versioned %>%
          filter(
            geo_value == this_loc, as_of == this_date,
            target_end_date <= (this_date - 2), !is.na(value)
          )
      } else {
        this_df <- flu_hosp_versioned %>%
          filter(
            geo_value == this_loc, as_of == this_date,
            target_end_date <= (this_date - 4), !is.na(value)
          )
      }
      if (nrow(this_df) < 4) next()
      meta <- list(
        location = name_map$location_number[which(name_map$geo_value == this_loc)],
        location_name = name_map$location_name[which(name_map$geo_value == this_loc)],
        forecast_date = as.Date(this_date, origin = "1970-01-01"),
        reference_date = as.Date(NA),
        target_end_date_base = as.Date(this_date, origin = "1970-01-01") + 5
      )
      q <- forecast_additive_one(this_df$value, meta, use_stl, use_log)
      if (!is.null(q)) out <- bind_rows(out, q)
    }
  }
  out
}

run_season_2324 <- function(use_stl, use_log) {
  out <- tibble()
  for (this_loc in unique(flu_hosp_versioned$geo_value)) {
    print(this_loc)
    for (this_date in forecast_dates_2324) {
      this_date <- as.Date(this_date, origin='1970-01-01')
      this_df <- flu_hosp_versioned %>%
        filter(geo_value == this_loc, as_of == this_date, target_end_date < this_date, !is.na(value))
      if (nrow(this_df) < 4) next()
      meta <- list(
        location = name_map$location_number[which(name_map$geo_value == this_loc)],
        location_name = name_map$location_name[which(name_map$geo_value == this_loc)],
        forecast_date = as.Date(this_date, origin = "1970-01-01"),
        reference_date = as.Date(this_date, origin = "1970-01-01") + 3,
        target_end_date_base = as.Date(this_date, origin = "1970-01-01") + 3
      )
      q <- forecast_additive_one(this_df$value, meta, use_stl, use_log)
      if (!is.null(q)) out <- bind_rows(out, q)
    }
  }
  out
}

run_season_2425 <- function(use_stl, use_log) {
  out <- tibble()
  for (this_loc in unique(df_hosp_2425$geo_value)) {
    print(this_loc)
    for (this_date in forecast_dates_2425) {
      this_date <- as.Date(this_date, origin='1970-01-01')
      this_df <- df_hosp_2425 %>%
        filter(geo_value == this_loc, as_of == this_date, target_end_date < this_date, !is.na(value))
      if (nrow(this_df) < 4) next()
      idx <- which(name_map$geo_value == this_loc)
      if (length(idx) == 0) next()
      meta <- list(
        location = name_map$location_number[idx],
        location_name = name_map$location_name[idx],
        forecast_date = as.Date(this_date, origin = "1970-01-01"),
        reference_date = as.Date(this_date, origin = "1970-01-01") + 3,
        target_end_date_base = as.Date(this_date, origin = "1970-01-01") + 3
      )
      q <- forecast_additive_one(this_df$value, meta, use_stl, use_log)
      if (!is.null(q)) out <- bind_rows(out, q)
    }
  }
  out
}


configs <- tribble(
  ~use_stl, ~use_log, ~out_file,
  FALSE, TRUE, "benchmark forecasts/trend baseline/naiveETS_damped_flu_hosp.parquet"
)
if (only_log) configs <- filter(configs, use_log, !use_stl)

for (i in seq_len(nrow(configs))) {
  use_stl <- configs$use_stl[[i]]
  use_log <- configs$use_log[[i]]
  path <- configs$out_file[[i]]
  message(
    "Additive AAN damped — STL=", use_stl, ", log(y+1)=", use_log,
    " -> ", path
  )
  combined <- bind_rows(
    run_season_2122(use_stl, use_log),
    run_season_2223(use_stl, use_log),
    run_season_2324(use_stl, use_log),
    run_season_2425(use_stl, use_log)
  )
  write_parquet(combined, path)
  message("  Wrote ", nrow(combined), " rows.")
  # keep the log / no-STL frame in memory for the reformatting stage below
  if (!use_stl && use_log) ets_log <- combined
}

message("Done (additive damped ETS).")


# ==========================================================================
# Reformatting: ETS output -> log-baseline schema (*_ets_models.parquet)
# Inlined from flu_weeklyhosp_format_conversion.R.
# Runs on the log / no-STL variant only -- the one the scoring pipeline uses.
# The unformatted parquets above are still written; this one is additional.
#
# Output schema: horizon, forecast_date, target_end_date, target,
#   abbreviation, quantile, value, location, location_name, Model
# Notes:
#   - forecast dates are NOT filtered to the logbaseline, so the Wednesday
#     2023-24 / 2024-25 forecasts (absent from the logbaseline) are kept;
#     the reference supplies only valid locations and the abbreviation map
#   - horizon parsed from target: "N wk ahead inc flu hosp" -> N (1-based)
#   - territories absent from the logbaseline (60/66/69) drop via valid_locs
# ==========================================================================

ETS_MODEL_LABEL <- "ETS Additive Damped - Log"
LBF <- "benchmark forecasts/random walk baseline/"
OUT <- "benchmark forecasts/trend baseline/"

if (!exists("ets_log")) {
  message("Reformatting skipped: log / no-STL config was not run.")
} else {

  ref <- read_parquet(paste0(LBF, "flu_hosp_log_baseline_h12.pq"))
  valid_locs <- unique(ref$location)
  abbrev_map <- ref %>% dplyr::select(location, abbreviation) %>% distinct()

  out_flu <- ets_log %>%
    filter(grepl("wk ahead inc flu hosp", target),
           location %in% valid_locs) %>%
    mutate(quantile        = round(quantile, 3),
           horizon         = as.integer(sub(" wk ahead.*", "", target)),
           forecast_date   = as.Date(forecast_date, origin='1970-01-01'),
           target_end_date = as.Date(target_end_date, origin='1970-01-01'),
           target          = "inc hosp",
           Model           = ETS_MODEL_LABEL) %>%
    left_join(abbrev_map, by = "location") %>%
    dplyr::select(horizon, forecast_date, target_end_date, target,
                  abbreviation, quantile, value, location, location_name, Model)
  write_parquet(out_flu, paste0(OUT, "flu_weeklyhosp_ets_models.parquet"))
  message("  flu_weeklyhosp_ets_models.parquet: ", nrow(out_flu), " rows.")

  message("Done (reformatted flu *_ets_models parquet).")
}
