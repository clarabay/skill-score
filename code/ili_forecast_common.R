# Shared ILI versioned data + FluSight 2015-2020 forecast dates.
# Used by ets_*_forecast_ili.R, tslm_forecast_ili.R, and naive ILI scripts.

require(tidyverse)
require(arrow)
require(lubridate)

#' @return list(ili_versioned, forecast_dates, seasons, ili_name_map, ilinet_all, targets)
build_ili_forecast_bundle <- function() {
  ilinet_all <- read_parquet("surveillance data/observed data/ilinet_truth.parquet")
  targets    <- read_csv("dat/targets_2015-2020.csv", show_col_types = FALSE)

  ili_versioned <- read_parquet("surveillance data/versioned_data/influenza_ili_new.parquet") %>%
    mutate(as_of = as.Date(as_of))

  forecast_dates <- read_csv("dat/ILI_ensemble_forecast_dates.csv",
                              show_col_types = FALSE) %>%
    mutate(
      ens_forecast_date = as.Date(ens_forecast_date),
      mmwr_year = epiyear(ens_forecast_date),
      mmwr_week = epiweek(ens_forecast_date),
      season    = ifelse(mmwr_week >= 40,
                         paste(mmwr_year, mmwr_year + 1, sep = "-"),
                         paste(mmwr_year - 1, mmwr_year, sep = "-"))
    ) %>%
    select(ens_forecast_date, season)

  seasons <- unique(targets$season)

  ili_name_map <- data.frame(
    abbreviation = c("nat", "hhs1", "hhs2", "hhs3", "hhs4", "hhs5",
                     "hhs6", "hhs7", "hhs8", "hhs9", "hhs10"),
    full = unique(ilinet_all$location)
  )

  list(
    ili_versioned  = ili_versioned,
    forecast_dates = forecast_dates,
    seasons        = seasons,
    ili_name_map   = ili_name_map,
    ilinet_all     = ilinet_all,
    targets        = targets
  )
}
