#!/usr/bin/env Rscript
# Additive damped ETS (AAN, damped=TRUE) for COVID deaths (weekly), cases (weekly), and hospitalizations
# (daily 1–28 day ahead for legacy + 2023–24; weekly 0–3 wk for 2024–25). Four outputs: (STL×log)×{on,off}.
# Training windows match naive_forecast_covid.R. Run from repository root.
#
# Performance: training series are pre-aggregated by (geo, forecast_date) to avoid repeated filter();
# results accumulated in lists then one bind_rows per stage.

rm(list = ls())

args      <- commandArgs(trailingOnly = TRUE)
only_log  <- "--only-log" %in% args

require(tidyverse)
require(arrow)
require(forecast)

source("code/create trend baseline/ets_forecast_covid_helpers.R")

name_map <- read.csv("dat/name_map.csv", stringsAsFactors = FALSE)
# Named vectors for O(1) location lookup
loc_num <- stats::setNames(name_map$location_number, name_map$geo_value)
loc_name <- stats::setNames(name_map$location_name, name_map$geo_value)

h_weekly   <- 12L
h_hosp2425 <- 12L
h_daily    <- 84L

# ── Deaths (weekly) ─────────────────────────────────────────────────────
df_death <- read_parquet("surveillance data/observed data/truth_deaths.parquet") %>%
  mutate(obs_value = ifelse(value < 0, NA, value)) %>%
  arrange(target_end_date)
forecast_dates_deaths <- unique(as.Date(df_death$target_end_date)) - 5

df_death_versioned <- read_parquet("surveillance data/versioned_data/covid_deaths.parquet")
df_death_versioned_agg <- df_death_versioned %>%
  mutate(value = ifelse(value < 0, NA, value)) %>%
  mutate(
    as_of = as.Date(as_of),
    target_end_date = ceiling_date(time_value, "week", week_start = 7) - 1
  ) %>%
  group_by(as_of, geo_value, target_end_date) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(geo_value)

death_ready <- df_death_versioned_agg %>%
  filter(target_end_date <= as_of - 2, !is.na(value)) %>%
  arrange(geo_value, as_of, target_end_date) %>%
  group_by(geo_value, as_of) %>%
  summarise(y_raw = list(value), .groups = "drop") %>%
  filter(lengths(y_raw) >= 4L, as_of %in% forecast_dates_deaths)

# ── Cases (weekly, states only) ───────────────────────────────────────────
df_case <- read_parquet("surveillance data/observed data/truth_cases.parquet") %>%
  mutate(obs_value = ifelse(value < 0, NA, value)) %>%
  filter(geo_type == "state") %>%
  arrange(target_end_date)
forecast_dates_cases <- unique(as.Date(df_case$target_end_date)) - 5

df_case_versioned <- read_parquet("surveillance data/versioned_data/covid_cases.parquet")
df_case_versioned_agg <- df_case_versioned %>%
  mutate(value = ifelse(value < 0, NA, value)) %>%
  mutate(
    as_of = as.Date(as_of),
    target_end_date = ceiling_date(time_value, "week", week_start = 7) - 1
  ) %>%
  group_by(as_of, geo_value, target_end_date) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(geo_value)

case_ready <- df_case_versioned_agg %>%
  filter(target_end_date <= as_of - 2, !is.na(value)) %>%
  arrange(geo_value, as_of, target_end_date) %>%
  group_by(geo_value, as_of) %>%
  summarise(y_raw = list(value), .groups = "drop") %>%
  filter(lengths(y_raw) >= 4L, as_of %in% forecast_dates_cases)

# ── Hospitalizations: daily (truth Mondays + 2023–24 window) ─────────────
df_hosp_truth <- read_parquet("surveillance data/observed data/truth_hosp.parquet") %>%
  mutate(obs_value = ifelse(value < 0, NA, value)) %>%
  arrange(target_end_date)
forecast_dates_hosp_mon <- unique(as.Date(df_hosp_truth$target_end_date[wday(df_hosp_truth$target_end_date) == 2]))

df_hosp_versioned <- read_parquet("surveillance data/versioned_data/covid_hosp.parquet") %>%
  mutate(as_of = as.Date(as_of), time_value = as.Date(time_value)) %>%
  arrange(geo_value)

df_hosp_2324_obs <- read_parquet('surveillance data/observed data/COVID19-observed.parquet') %>%
  filter(outcome == 'hosp')  %>%
  mutate(
    target_end_date = as.Date(target_end_date),
    obs_value = ifelse(obs_value < 0, NA, obs_value)) %>%
  arrange(target_end_date)

forecast_dates_2324 <- unique(df_hosp_2324_obs$target_end_date[wday(df_hosp_2324_obs$target_end_date) == 2])
forecast_dates_2324 <- forecast_dates_2324[
  as.Date("2023-01-01") < forecast_dates_2324 & forecast_dates_2324 < as.Date("2024-05-01")
]

hosp_daily_dates <- unique(c(forecast_dates_hosp_mon, forecast_dates_2324))
hosp_daily_ready <- df_hosp_versioned %>%
  filter(time_value < as_of, !is.na(value)) %>%
  arrange(geo_value, as_of, time_value) %>%
  group_by(geo_value, as_of) %>%
  summarise(y_raw = list(value), .groups = "drop") %>%
  filter(lengths(y_raw) >= 4L, as_of %in% hosp_daily_dates)

# ── Hospitalizations: 2024–25 weekly (Wednesday issue dates) ──────────────
df_hosp_2425 <- read_parquet("surveillance data/versioned_data/covid_hosp_24-25.parquet") %>%
  mutate(
    obs_value = ifelse(value < 0, NA, value),
    as_of = as.Date(as_of),
    target_end_date = as.Date(target_end_date)
  ) %>%
  filter(target_end_date <= as.Date("2024-05-04") | as.Date("2024-11-09") <= target_end_date) %>%
  arrange(target_end_date)

forecast_dates_2425 <- unique(df_hosp_2425$as_of)
forecast_dates_2425 <- forecast_dates_2425[as.Date("2024-11-20") <= forecast_dates_2425]

hosp2425_ready <- df_hosp_2425 %>%
  filter(target_end_date < as_of, !is.na(value)) %>%
  arrange(geo_value, as_of, target_end_date) %>%
  group_by(geo_value, as_of) %>%
  summarise(y_raw = list(value), .groups = "drop") %>%
  filter(lengths(y_raw) >= 4L, as_of %in% forecast_dates_2425)

meta_or_skip <- function(geo) {
  if (is.null(loc_num[[geo]])) return(NULL)
  list(
    location = loc_num[[geo]],
    location_name = loc_name[[geo]]
  )
}

run_deaths <- function(use_stl, use_log) {
  chunks <- list()
  k <- 0L
  nr <- nrow(death_ready)
  message("    deaths: ", nr, " location×dates")
  for (ri in seq_len(nr)) {
    if (ri %% 400L == 1L) message("      deaths ", ri, "/", nr)
    row <- death_ready[ri, ]
    y_raw <- row$y_raw[[1]]
    this_loc <- row$geo_value
    this_week <- row$as_of
    nm <- meta_or_skip(this_loc)
    if (is.null(nm)) next()
    meta <- c(
      nm,
      list(
        forecast_date = this_week,
        reference_date = as.Date(NA),
        outcome = "death",
        target_end_date_base = this_week + 5
      )
    )
    pred_fn <- function(pred, meta, ls, adj) {
      pred_to_quant_covid_weekly(pred, meta, ls, adj, h_weekly, "death")
    }
    q <- forecast_additive_covid_one(
      y_raw, meta, use_stl, use_log, h_weekly, STL_PERIOD_WEEKLY, pred_fn
    )
    if (!is.null(q)) {
      k <- k + 1L
      chunks[[k]] <- q
    }
  }
  if (k == 0L) tibble() else bind_rows(chunks)
}

run_cases <- function(use_stl, use_log) {
  chunks <- list()
  k <- 0L
  nr <- nrow(case_ready)
  message("    cases: ", nr, " location×dates")
  for (ri in seq_len(nr)) {
    if (ri %% 400L == 1L) message("      cases ", ri, "/", nr)
    row <- case_ready[ri, ]
    y_raw <- row$y_raw[[1]]
    this_loc <- row$geo_value
    this_week <- row$as_of
    nm <- meta_or_skip(this_loc)
    if (is.null(nm)) next()
    meta <- c(
      nm,
      list(
        forecast_date = this_week,
        reference_date = as.Date(NA),
        outcome = "case",
        target_end_date_base = this_week + 5
      )
    )
    pred_fn <- function(pred, meta, ls, adj) {
      pred_to_quant_covid_weekly(pred, meta, ls, adj, h_weekly, "cases")
    }
    q <- forecast_additive_covid_one(
      y_raw, meta, use_stl, use_log, h_weekly, STL_PERIOD_WEEKLY, pred_fn
    )
    if (!is.null(q)) {
      k <- k + 1L
      chunks[[k]] <- q
    }
  }
  if (k == 0L) tibble() else bind_rows(chunks)
}

run_hosp_daily <- function(use_stl, use_log) {
  chunks <- list()
  k <- 0L
  nr <- nrow(hosp_daily_ready)
  message("    hosp daily: ", nr, " location×dates")
  for (ri in seq_len(nr)) {
    if (ri %% 500L == 1L) message("      hosp daily ", ri, "/", nr)
    row <- hosp_daily_ready[ri, ]
    y_raw <- row$y_raw[[1]]
    this_loc <- row$geo_value
    this_date <- row$as_of
    nm <- meta_or_skip(this_loc)
    if (is.null(nm)) next()
    meta <- c(
      nm,
      list(
        forecast_date = this_date,
        reference_date = as.Date(NA),
        outcome = "hosp"
      )
    )
    pred_fn <- function(pred, meta, ls, adj) {
      pred_to_quant_covid_daily_hosp(pred, meta, ls, adj, h_daily)
    }
    q <- forecast_additive_covid_one(
      y_raw, meta, use_stl, use_log, h_daily, STL_PERIOD_DAILY, pred_fn
    )
    if (!is.null(q)) {
      k <- k + 1L
      chunks[[k]] <- q
    }
  }
  if (k == 0L) tibble() else bind_rows(chunks)
}

run_hosp_2425 <- function(use_stl, use_log) {
  chunks <- list()
  k <- 0L
  nr <- nrow(hosp2425_ready)
  message("    hosp 2024–25: ", nr, " location×dates")
  for (ri in seq_len(nr)) {
    if (ri %% 200L == 1L) message("      hosp 2425 ", ri, "/", nr)
    row <- hosp2425_ready[ri, ]
    y_raw <- row$y_raw[[1]]
    this_loc <- row$geo_value
    this_date <- row$as_of
    nm <- meta_or_skip(this_loc)
    if (is.null(nm)) next()
    meta <- c(
      nm,
      list(
        forecast_date = this_date,
        reference_date = this_date + 3,
        outcome = "hosp"
      )
    )
    pred_fn <- function(pred, meta, ls, adj) {
      pred_to_quant_covid_hosp2425(pred, meta, ls, adj, h_hosp2425)
    }
    q <- forecast_additive_covid_one(
      y_raw, meta, use_stl, use_log, h_hosp2425, STL_PERIOD_WEEKLY, pred_fn
    )
    if (!is.null(q)) {
      k <- k + 1L
      chunks[[k]] <- q
    }
  }
  if (k == 0L) tibble() else bind_rows(chunks)
}


configs <- tribble(
  ~use_stl, ~use_log, ~out_file,
  FALSE, TRUE, "benchmark forecasts/trend baseline/ets_additive_damped_covid_log.parquet",
 )

if (only_log) configs <- filter(configs, use_log, !use_stl)

for (i in seq_len(nrow(configs))) {
  use_stl <- configs$use_stl[[i]]
  use_log <- configs$use_log[[i]]
  path <- configs$out_file[[i]]
  message("COVID additive ETS — STL=", use_stl, ", log(y+1)=", use_log, " -> ", path)
  combined <- bind_rows(
    run_deaths(use_stl, use_log),
    run_cases(use_stl, use_log),
    run_hosp_daily(use_stl, use_log),
    run_hosp_2425(use_stl, use_log)
  )
  write_parquet(combined, path)
  message("  Wrote ", nrow(combined), " rows.")
  # keep the log / no-STL frame in memory for the reformatting stage below
  if (!use_stl && use_log) ets_log <- combined
}

message("Done (additive damped ETS COVID).")


# ==========================================================================
# Reformatting: ETS output  (*_ets_models.parquet)
#
# Output schema: horizon, forecast_date|reference_date, target_end_date,
#   target, abbreviation, quantile, value, location, location_name, Model
# Territories absent from the logbaseline (60/66/69/78) drop out via valid_locs.
# ==========================================================================

ETS_MODEL_LABEL <- "ETS Additive Damped - Log"
LBF <- "benchmark forecasts/random walk baseline/"
OUT <- "benchmark forecasts/trend baseline/"

if (!exists("ets_log")) {
  message("Reformatting skipped: log / no-STL config was not run.")
} else {

  # ---- COVID cases: horizon parsed from target string ---------------------
  ref <- read_parquet(paste0(LBF, "covid_case_log_baseline_h12.pq"))
  valid_dates  <- as.Date(unique(ref$forecast_date))
  valid_locs   <- unique(ref$location)
  valid_quants <- unique(ref$quantile)
  abbrev_map   <- ref %>% dplyr::select(location, abbreviation) %>% distinct()

  out_case <- ets_log %>%
    filter(outcome == "case",
           location %in% valid_locs,
           quantile %in% valid_quants,
           as.Date(forecast_date) %in% valid_dates) %>%
    mutate(quantile      = round(quantile, 3),
           horizon       = as.integer(sub(" wk ahead.*", "", target)),
           forecast_date = as.Date(forecast_date),
           target        = "inc case",
           Model         = ETS_MODEL_LABEL) %>%
    left_join(abbrev_map, by = "location") %>%
    dplyr::select(horizon, forecast_date, target_end_date, target,
                  abbreviation, quantile, value, location, location_name, Model)
  write_parquet(out_case, paste0(OUT, "covid_case_ets_models.parquet"))
  message("  covid_case_ets_models.parquet: ", nrow(out_case), " rows.")

  # ---- COVID deaths: --------
  ref <- read_parquet(paste0(LBF, "covid_death_log_baseline_h12.pq"))
  valid_locs  <- unique(ref$location)
  abbrev_map  <- ref %>% dplyr::select(location, abbreviation) %>% distinct()
  horizon_map <- ref %>%
    dplyr::select(forecast_date, target_end_date, horizon) %>%
    mutate(forecast_date   = as.Date(forecast_date),
           target_end_date = as.Date(target_end_date)) %>%
    distinct()

  out_death <- ets_log %>%
    filter(outcome == "death", location %in% valid_locs) %>%
    mutate(quantile        = round(quantile, 3),
           forecast_date   = as.Date(forecast_date),
           target_end_date = as.Date(target_end_date),
           target          = "inc death",
           Model           = ETS_MODEL_LABEL) %>%
    inner_join(horizon_map, by = c("forecast_date", "target_end_date")) %>%
    left_join(abbrev_map, by = "location") %>%
    dplyr::select(horizon, forecast_date, target_end_date, target,
                  abbreviation, quantile, value, location, location_name, Model)
  write_parquet(out_death, paste0(OUT, "covid_death_ets_models.parquet"))
  message("  covid_death_ets_models.parquet: ", nrow(out_death), " rows.")

  # ---- COVID daily hosp: inner_join drops ETS "1 day ahead" (Tuesday) ----
  ref <- read_parquet(paste0(LBF, "covid_dailyhosp_log_baseline_h12.pq"))
  valid_dates <- as.Date(unique(ref$forecast_date))
  valid_locs  <- unique(ref$location)
  abbrev_map  <- ref %>% dplyr::select(location, abbreviation) %>% distinct()
  horizon_map <- ref %>%
    dplyr::select(forecast_date, target_end_date, horizon) %>%
    mutate(forecast_date   = as.Date(forecast_date),
           target_end_date = as.Date(target_end_date)) %>%
    distinct()

  out_dailyhosp <- ets_log %>%
    filter(outcome == "hosp",
           grepl("day ahead", target),
           location %in% valid_locs,
           as.Date(forecast_date) %in% valid_dates) %>%
    mutate(quantile        = round(quantile, 3),
           forecast_date   = as.Date(forecast_date),
           target_end_date = as.Date(target_end_date),
           target          = "inc hosp",
           Model           = ETS_MODEL_LABEL) %>%
    inner_join(horizon_map, by = c("forecast_date", "target_end_date")) %>%
    left_join(abbrev_map, by = "location") %>%
    dplyr::select(horizon, forecast_date, target_end_date, target,
                  abbreviation, quantile, value, location, location_name, Model)
  write_parquet(out_dailyhosp, paste0(OUT, "covid_dailyhosp_ets_models.parquet"))
  message("  covid_dailyhosp_ets_models.parquet: ", nrow(out_dailyhosp), " rows.")

  # ---- COVID weekly hosp: join on (reference_date, ted); no date filter --
  
  ref <- read_parquet(paste0(LBF, "covid_hosp_log_baseline_2025_h12.pq"))
  valid_locs  <- unique(ref$location)
  abbrev_map  <- ref %>% dplyr::select(location, abbreviation) %>% distinct()
  horizon_map <- ref %>%
    dplyr::select(reference_date, target_end_date, horizon) %>%
    mutate(reference_date  = as.Date(reference_date),
           target_end_date = as.Date(target_end_date)) %>%
    distinct()

  out_weeklyhosp <- ets_log %>%
    filter(outcome == "hosp",
           grepl("wk ahead", target),
           location %in% valid_locs) %>%
    mutate(quantile        = round(quantile, 3),
           reference_date  = as.Date(reference_date),
           target_end_date = as.Date(target_end_date),
           target          = "inc hosp",
           Model           = ETS_MODEL_LABEL) %>%
    inner_join(horizon_map, by = c("reference_date", "target_end_date")) %>%
    left_join(abbrev_map, by = "location") %>%
    dplyr::select(horizon, reference_date, target_end_date, target,
                  abbreviation, quantile, value, location, location_name, Model)
  write_parquet(out_weeklyhosp, paste0(OUT, "covid_weeklyhosp_ets_models.parquet"))
  message("  covid_weeklyhosp_ets_models.parquet: ", nrow(out_weeklyhosp), " rows.")

  message("Done (reformatted COVID *_ets_models parquets).")
}
