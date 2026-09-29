rm(list=ls())

require(tidyverse)
require(purrr)
require(arrow)
require(MMWRweek)

source('code/forecast_functions.R')

ilinet_all <- read_parquet('surveillance data/observed data/ilinet_truth.parquet')

season_targets <- c("1 wk ahead", "2 wk ahead", "3 wk ahead", "4 wk ahead")

targets <- read_csv('dat/targets_2015-2020.csv')

forecast_dates <- bind_rows(
  tibble(
    year = 2015,
    epi_week = c(1:(21 + 4), 42:52)),
  tibble(
    year = 2016,
    epi_week = c(1:(18 + 4), 43:52)),
  tibble(
    year = 2017,
    epi_week = c(1:(18 + 4), 43:52)),
  tibble(
    year = 2018,
    epi_week = c(1:(18 + 4), 42:52)),
  tibble(
    year = 2019,
    epi_week = c(1:(18 + 4), 42:52)),
  tibble(
    year = 2020,
    epi_week = 1:(9 + 4))) %>%
  mutate(
    season = ifelse(epi_week >= 40, paste(year, year+1, sep='-'), paste(year-1, year, sep='-')),
    epi_week_start_date = MMWRweek2Date(year, epi_week, MMWRday=1),
    target_end_date = epi_week_start_date + 6)

quantiles23 <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
  0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

seasons <- unique(targets$season)

hindcasts_param <- tibble()
hindcasts_bin <- tibble()
hindcasts_quant <- tibble()
for (this_loc in unique(ilinet_all$location)) {
#  for (this_target in season_targets) {
    for (this_season in seasons) {
#      this_loc <- 'US National'; this_target <- '1 wk ahead'; this_season = '2016-2017'; this_week_start = as.Date('2016-10-23')
      season_weeks <- filter(forecast_dates, season == this_season)
      for (this_week_start in season_weeks$epi_week_start_date) {
        this_week_start <- as.Date(this_week_start, origin='1970-01-01')
        this_df <- filter(ilinet_all, location == this_loc, 
          (this_week_start - 7) <= week_start, week_start <= (this_week_start + 7))
        #      beta_fit <- fit_beta(this_df$weighted_ili, range=c(0, 100))
        logitnormal_fit <- fit_logitnorm(this_df$weighted_ili/100, weights=c(1, 2, 1))
        #      plot(density(rbeta(1e4, shape1=beta_fit[['shape1']], shape2=beta_fit[['shape2']])))
        #      lines(density(boot::inv.logit(rnorm(1e4, logitnormal_fit[[1]], logitnormal_fit[[2]]))), col='red')
        this_param <- tibble(
          location = this_loc,
          target = 'wk ahead', #this_target,
          season = this_season,
          target_end_date = as.Date(this_week_start + 6, origin='1970-01-01'),
          logitnorm_mean = logitnormal_fit$logit_mean,
          logitnorm_sd = logitnormal_fit$logit_sd
        )
        these_bins <- filter(targets, season == this_season, target == '1 wk ahead', 
          location == this_loc, !is.na(bin_end_notincl)) %>%
          mutate(
            bin_start_incl = as.numeric(bin_start_incl)/100,
            bin_end_notincl = round(as.numeric(bin_end_notincl), 1)/100
          ) %>%
          unique()
        if (max(these_bins$bin_end_notincl) != 1) stop("max bin end value should be 1")
        this_bin <- this_param %>%
          mutate(bv = map2(logitnorm_mean, logitnorm_sd, ~ logitnorm_to_bin(logit_mean=.x, logit_sd=.y, 
            bins=these_bins))) %>%
          unnest(bv) %>%
          dplyr::select(-logitnorm_mean, -logitnorm_sd)
        if (!all.equal(sum(this_bin$value), 1)) stop('bin values did not add to one')
        this_quant <- this_param %>%
          mutate(bv = map2(logitnorm_mean, logitnorm_sd, ~ logitnorm_to_quant(logit_mean=.x, logit_sd=.y, 
            quantiles=quantiles23))) %>%
          unnest(bv) %>%
          dplyr::select(-logitnorm_mean, -logitnorm_sd)
        hindcasts_param <- bind_rows(hindcasts_param, this_param)
        hindcasts_bin <- bind_rows(hindcasts_bin, this_bin)
        hindcasts_quant <- bind_rows(hindcasts_quant, this_quant)
      }
    }
#  }
}

write_parquet(hindcasts_bin, 'benchmark forecasts/hindcasts/flusight-hindcast-Bin.parquet')
