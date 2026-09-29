rm(list=ls())

require(tidyverse)
require(purrr)
require(arrow)
require(MMWRweek)

source('code/forecast_functions.R')

horizons_to_forecast <- 12

quantiles <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 
  0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)


# ILI
source("code/ili_forecast_common.R")
bundle         <- build_ili_forecast_bundle()
ili_versioned  <- bundle$ili_versioned
forecast_dates <- bundle$forecast_dates
seasons        <- bundle$seasons
ili_name_map   <- bundle$ili_name_map
ilinet_all     <- bundle$ilinet_all
targets        <- bundle$targets

quantiles23 <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
  0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

short_term_targets <- c("1 wk ahead", "2 wk ahead", "3 wk ahead", "4 wk ahead")
if (horizons_to_forecast != 4) {
  short_term_targets <- paste0(1:horizons_to_forecast, " wk ahead")
}

naive_ili_fcasts_bin <- tibble()
naive_ili_fcasts_quant <- tibble()
for (this_loc in unique(ilinet_all$location)) {
  print(this_loc)
  loc_abbr <- ili_name_map$abbreviation[which(ili_name_map$full == this_loc)]
  for (i_date in seq_len(nrow(forecast_dates))) {
    this_forecast_date <- forecast_dates$ens_forecast_date[i_date]
    this_season        <- forecast_dates$season[i_date]
    these_bins <- filter(targets, season == this_season, target == short_term_targets[1],
      location == this_loc, !is.na(bin_end_notincl)) %>%
      mutate(
        bin_start_incl  = as.numeric(bin_start_incl) / 100,
        bin_end_notincl = round(as.numeric(bin_end_notincl), 1) / 100
      ) %>%
      unique()
    if (max(these_bins$bin_end_notincl) != 1) stop("max bin end value should be 1")
    this_df <- filter(ili_versioned, geo_value == loc_abbr,
                      as_of == this_forecast_date, !is.na(wili)) %>%
      arrange(target_end_date)
    if (nrow(this_df) < 3) next()
    logitnormal_fit <- fit_logitnorm(this_df$wili / 100, zero_val = 1/1000)
    epi_week_sunday <- floor_date(this_forecast_date, "week", week_start = 7)
    this_param <- tibble(
        location        = this_loc,
        forecast_date   = this_forecast_date,
        target          = short_term_targets,
        target_end_date = as.Date(epi_week_sunday + 7 * ((1:horizons_to_forecast) - 1) - 1),
        logitnorm_mean  = logitnormal_fit$logit_mean,
        logitnorm_sd    = logitnormal_fit$logit_sd
      )
    this_bin <- this_param %>%
      mutate(bv = map2(logitnorm_mean, logitnorm_sd,
        ~ logitnorm_to_bin(logit_mean = .x, logit_sd = .y, bins = these_bins))) %>%
      unnest(bv) %>%
      dplyr::select(-logitnorm_mean, -logitnorm_sd)
    check_sums <- group_by(this_bin, target_end_date) %>%
      summarize(sum = sum(value)) %>% pull(sum)
    this_quant <- this_param %>%
      mutate(bv = map2(logitnorm_mean, logitnorm_sd,
        ~ logitnorm_to_quant(logit_mean = .x, logit_sd = .y, quantiles = quantiles23))) %>%
      unnest(bv) %>%
      dplyr::select(-logitnorm_mean, -logitnorm_sd)
    if (!all.equal(check_sums, rep(1, horizons_to_forecast))) stop('bin values did not add to one')
    naive_ili_fcasts_bin   <- bind_rows(naive_ili_fcasts_bin,   this_bin)
    naive_ili_fcasts_quant <- bind_rows(naive_ili_fcasts_quant, this_quant)
  }
}

if (horizons_to_forecast == 4) {
  write_parquet(naive_ili_fcasts_bin, 'benchmark forecasts/naive/naive_flu_ili.parquet')
} else {
  write_parquet(naive_ili_fcasts_bin, paste0('benchmark forecasts/naive/naive_flu_ili_h',
    horizons_to_forecast, '_2014-20.parquet'))
}


