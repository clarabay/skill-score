#!/usr/bin/env Rscript
# Additive damped ETS (AAN, damped=TRUE) for ILI using logit(weighted_ili/100).
# Training history is restricted to 2003/2004 season onward.
# Output: non-STL logit forecast parquet.

rm(list = ls())

require(tidyverse)
require(arrow)
require(forecast)
require(boot)

source("code/create trend baseline/ets_forecast_ili_logit_helpers.R")
source("code/forecast_functions.R")

bundle <- build_ili_forecast_bundle()
ili_versioned      <- bundle$ili_versioned
forecast_dates_tbl <- bundle$forecast_dates
ili_name_map       <- bundle$ili_name_map
ilinet_locations   <- unique(bundle$ilinet_all$location)
targets            <- bundle$targets

# 2003/2004 season starts at MMWR week 40 of 2003.
min_week_start <- as.Date("2003-09-28")

forecast_additive_logit_one <- function(y_raw, meta, bins = NULL) {
  if (length(y_raw) < 4L) return(NULL)

  x <- boot::logit(ili_to_prop(y_raw))
  if (any(!is.finite(x))) return(NULL)

  m <- tryCatch(ets(x, model = "AAN", damped = TRUE), error = function(e) NULL)
  if (is.null(m)) return(NULL)

  pred <- tryCatch(
    forecast(m, h = horizons_to_forecast, level = forecast_levels),
    error = function(e) NULL
  )
  if (is.null(pred)) return(NULL)

  q <- pred_to_quant_logit(pred, meta, forecast_adjustment = rep(0, horizons_to_forecast))
  b <- if (!is.null(bins))
    pred_to_bin_logit(pred, meta, bins, forecast_adjustment = rep(0, horizons_to_forecast))
  else NULL
  list(quant = q, bins = b)
}

run_all_ili_logit_2003on <- function() {
  out_quant <- tibble()
  out_bins  <- tibble()
  for (this_loc in ilinet_locations) {
    print(this_loc)
    loc_abbr <- ili_name_map$abbreviation[which(ili_name_map$full == this_loc)]
    for (this_forecast_date in forecast_dates_tbl$ens_forecast_date) {
      this_forecast_date <- as.Date(this_forecast_date, origin='1970-01-01')
      this_df <- filter(ili_versioned, geo_value == loc_abbr,
                        as_of == this_forecast_date,
                        target_end_date >= min_week_start,
                        !is.na(wili)) %>%
        arrange(target_end_date) %>%
        rename(weighted_ili = wili, week_start = target_end_date)

      y_raw <- this_df$weighted_ili
      if (length(y_raw) < 4L) next()

      epi_week_sunday <- floor_date(this_forecast_date, "week", week_start = 7)
      meta <- list(
        location             = this_loc,
        location_name        = this_loc,
        forecast_date        = this_forecast_date,
        reference_date       = as.Date(NA),
        target_end_date_base = epi_week_sunday - 1L
      )

      this_season <- forecast_dates_tbl$season[forecast_dates_tbl$ens_forecast_date == this_forecast_date]
      these_bins <- filter(targets,
                           season == this_season,
                           target == "1 wk ahead",
                           location == this_loc,
                           !is.na(bin_end_notincl)) %>%
        mutate(
          bin_start_incl  = as.numeric(bin_start_incl) / 100,
          bin_end_notincl = round(as.numeric(bin_end_notincl), 1) / 100
        ) %>%
        unique()

      result <- forecast_additive_logit_one(y_raw, meta,
                                            bins = if (nrow(these_bins) > 0) these_bins else NULL)
      if (!is.null(result)) {
        if (!is.null(result$quant)) out_quant <- bind_rows(out_quant, result$quant)
        if (!is.null(result$bins))  out_bins  <- bind_rows(out_bins,  result$bins)
      }
    }
  }
  list(quantiles = out_quant, bins = out_bins)
}

out_file_bins <- "benchmark forecasts/trend baseline/ets_additive_damped_ili_logit_2003on_bins_2014-20.parquet"
#message("Additive AAN damped (ILI) - logit(wILI/100), training >= 2003-09-28 -> ", out_file)
combined <- run_all_ili_logit_2003on()
#write_parquet(combined$quantiles, out_file)
#message("Wrote ", nrow(combined$quantiles), " rows to quantiles file.")
write_parquet(combined$bins, out_file_bins)
message("Wrote ", nrow(combined$bins), " rows to bins file.")
