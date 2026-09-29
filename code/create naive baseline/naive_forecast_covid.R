rm(list=ls())

require(tidyverse)
require(zoo)
require(fitdistrplus)
require(purrr)
require(arrow)

source('code/forecast_functions.R')
name_map <- read.csv("dat/name_map.csv")

horizons_to_forecast <- 12 # weeks

### deaths
df_death <- read_parquet('surveillance data/observed data/truth_deaths.parquet') %>%
  mutate(obs_value = ifelse(value < 0, NA, value)) %>%
  arrange(target_end_date)
forecast_dates_deaths <- unique(df_death$target_end_date) - 5 

quantiles_d <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
                 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

library(foreach)
library(doParallel)
library(lubridate) 

# Set up parallel backend
num_cores <- detectCores() - 1  # Use one less than available cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Make sure to export all required functions and objects to the workers
packages <- c("dplyr", "lubridate", "tibble", "fitdistrplus", 'tidyr', 'purrr')

# Load packages on all worker nodes
clusterEvalQ(cl, {
  library(dplyr)
  library(tibble)
  library(fitdistrplus)
  library(lubridate)
  library(purrr)
  library(tidyr)
})

df_death_versioned <- read_parquet('surveillance data/versioned_data/covid_deaths.parquet')
df_death_versioned_agg <- df_death_versioned %>%
  mutate(value = ifelse(value < 0, NA, value)) %>%
  mutate(target_end_date = ceiling_date(time_value, "week", week_start = 7)) %>%
  mutate(target_end_date = target_end_date-1) %>%
  group_by(as_of, geo_value, target_end_date) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = 'drop') %>%
  arrange(geo_value)

# Parallelize the outer loop
results_death <- foreach(this_loc = unique(df_death_versioned$geo_value), 
           .packages = c("dplyr", "lubridate", "tibble", "fitdistrplus", 'tidyr', 'purrr'),
           .combine = "bind_rows") %dopar% {
             
             # For progress tracking
             message(cat("Processing location:", this_loc, "\n"))
             
             loc_quant <- tibble()
             
             for (this_week in forecast_dates_deaths) {
               this_week <- as.Date(this_week, origin = "1970-01-01")
               this_df <- df_death_versioned_agg %>%
                 filter(as_of == this_week, geo_value == this_loc ,
                   target_end_date <= (this_week - 2), !is.na(value))
               if (nrow(this_df) < 3) next()
               nb_fit <- fit_nbinom(this_df$value, if_all_zeros='add1')
               this_param <- tibble(
                 outcome = 'death',
                 location = name_map$location_number[which(name_map$geo_value == this_loc)],
                 location_name = name_map$location_name[which(name_map$geo_value == this_loc)],
                 forecast_date = as.Date(this_week, origin='1970-01-01'),
                 horizon = (1:horizons_to_forecast) - 1,
                 target = paste0(1:horizons_to_forecast, ' wk ahead inc death'),
                 target_end_date = as.Date(this_week + 5, origin='1970-01-01') + ((1:horizons_to_forecast) - 1)*7,
                 nbinom_mu = nb_fit$mu,
                 nbinom_size = nb_fit$size
               )
               this_quant <- this_param %>%
                 mutate(qv = map2(nbinom_mu, nbinom_size, ~ nbinom_to_quant(mu=.x, size=.y, quantiles=quantiles_d))) %>%
                 unnest(qv) %>%
                 dplyr::select(-nbinom_mu, -nbinom_size)
               loc_quant <- bind_rows(loc_quant, this_quant)
             }
             
             # Return results for this location
             loc_quant
           }

stopCluster(cl)
naive_d_quant <- results_death

if (horizons_to_forecast == 4) {
  write_parquet(naive_d_quant, 'benchmark forecasts/naive/naive_covid_deaths.parquet')
} else {
  write_parquet(naive_d_quant, 
    paste0('benchmark forecasts/naive/naive_covid_deaths_h', horizons_to_forecast, '.parquet'))
}

################################################################
### cases
df_case <- read_parquet('surveillance data/observed data/truth_cases.parquet') %>%
  mutate(obs_value = ifelse(value < 0, NA, value)) %>%
  arrange(target_end_date) %>%
  filter(geo_type == 'state')
forecast_dates_cases <- unique(df_case$target_end_date) -5 

quantiles_c <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)

df_case_versioned <- read_parquet('surveillance data/versioned_data/covid_cases.parquet')
df_case_versioned_agg <- df_case_versioned %>%
  mutate(value = ifelse(value < 0, NA, value)) %>%
  mutate(target_end_date = ceiling_date(time_value, "week", week_start = 7)) %>%
  mutate(target_end_date = target_end_date-1) %>%
  group_by(as_of, geo_value, target_end_date) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = 'drop') %>%
  arrange(geo_value)

# Set up parallel backend
num_cores <- detectCores() - 1  # Use one less than available cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Make sure to export all required functions and objects to the workers
packages <- c("dplyr", "lubridate", "tibble", "fitdistrplus", 'tidyr', 'purrr')

# Load packages on all worker nodes
clusterEvalQ(cl, {
  library(dplyr)
  library(tibble)
  library(fitdistrplus)
  library(lubridate)
  library(purrr)
  library(tidyr)
})
# Parallelize the outer loop
results_cases <- foreach(this_loc = unique(df_case_versioned$geo_value), 
    .packages = c("dplyr", "lubridate", "tibble", "fitdistrplus", 'tidyr', 'purrr'),
    .combine = "bind_rows") %dopar% {
     
     # For progress tracking
     message(cat("Processing location:", this_loc, "\n"))
     
     loc_quant <- tibble()
     
     for (this_week in forecast_dates_cases) {
       this_week <- as.Date(this_week, origin = "1970-01-01")
       
       this_df <- df_case_versioned_agg %>%
         filter(as_of == this_week, geo_value == this_loc ,
                target_end_date <= (this_week -2), !is.na(value)) %>%
         mutate(lobs = log(value + 1))
       if (nrow(this_df) < 3) next()
       nb_fit <- fit_nbinom(this_df$value, if_all_zeros='add1')
       this_quant <- tibble(
         outcome = 'case',
         location = name_map$location_number[which(name_map$geo_value == this_loc)],
         location_name = name_map$location_name[which(name_map$geo_value == this_loc)],
         forecast_date = as.Date(this_week, origin='1970-01-01'),
         horizon = (1:horizons_to_forecast) - 1,
         target = paste0(1:horizons_to_forecast, ' wk ahead inc death'),
         target_end_date = as.Date(this_week + 5, origin='1970-01-01') + ((1:horizons_to_forecast) - 1)*7,
         nbinom_mu = nb_fit$mu,
         nbinom_size = nb_fit$size
       ) %>%
         mutate(qv = map2(nbinom_mu, nbinom_size, ~ nbinom_to_quant(mu=.x, size=.y, quantiles=quantiles_c))) %>%
         unnest(qv) %>%
         dplyr::select(-nbinom_mu, -nbinom_size)
       loc_quant <- bind_rows(loc_quant, this_quant)
     }
     
     # Return results for this location
     loc_quant
    }

# Stop the cluster when done
stopCluster(cl)

# Results are now in the 'results_death' dataframe
naive_c_quant <- results_cases

if (horizons_to_forecast == 4) {
  write_parquet(naive_c_quant, 'benchmark forecasts/naive/naive_covid_cases.parquet')
} else {
  write_parquet(naive_c_quant, 
    paste0('benchmark forecasts/naive/naive_covid_cases_h', horizons_to_forecast, '.parquet'))
}

###########################################
### hospitalizations
df_hosp <- read_parquet('surveillance data/observed data/truth_hosp.parquet') %>%
  mutate(obs_value = ifelse(value < 0, NA, value)) %>%
  arrange(target_end_date)

quantiles_h <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
                 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

# select Mondays
forecast_dates <- unique(df_hosp$target_end_date[wday(df_hosp$target_end_date) == '2'])

df_hosp_versioned <- read_parquet('surveillance data/versioned_data/covid_hosp.parquet') %>% 
  arrange(geo_value)

if (horizons_to_forecast == 4) {
  these_target_days <- 1:28
} else {
  these_target_days <- (1:horizons_to_forecast) * 7 - 2
}

naive_h_quant_list <- list()  # Use a list to collect tibbles
i <- 1  # Counter for list indexing

for (this_loc in unique(df_hosp_versioned$geo_value)) {
  print(this_loc)
  for (this_date in forecast_dates) {
    this_date <- as.Date(this_date, origin = "1970-01-01")
    this_df <- df_hosp_versioned %>%
      filter(as_of == this_date, geo_value == this_loc,
             time_value < this_date, !is.na(value))
    if (nrow(this_df) < 3) next()
    nb_fit <- fit_nbinom(this_df$value, if_all_zeros='add1')
    
    this_quant <- tibble(
      outcome = 'hosp',
      location = name_map$location_number[which(name_map$geo_value == this_loc)],
      location_name = name_map$location_name[which(name_map$geo_value == this_loc)],
      forecast_date = this_date,
      horizon = these_target_days/7 - 1,
      target = paste0(these_target_days, ' day ahead inc hosp'),
      target_end_date = this_date + these_target_days,
      nbinom_mu = nb_fit$mu,
      nbinom_size = nb_fit$size
    ) %>%
      mutate(qv = map2(nbinom_mu, nbinom_size, ~ nbinom_to_quant(mu=.x, size=.y, quantiles=quantiles_h))) %>%
      unnest(qv) %>%
      dplyr::select(-nbinom_mu, -nbinom_size)
    
    naive_h_quant_list[[i]] <- this_quant
    i <- i + 1
  }
}

naive_h_quant <- bind_rows(naive_h_quant_list)

###########################################
### hospitalizations 2023-2024
df_hosp_2324 <- read_parquet('surveillance data/observed data/COVID19-observed.parquet') %>%
  filter(outcome == 'hosp')  %>%
  mutate(
    target_end_date = as.Date(target_end_date),
    obs_value = ifelse(obs_value < 0, NA, obs_value)) %>%
  arrange(target_end_date)

quantiles_h <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
                 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

# select Mondays
forecast_dates <- unique(df_hosp_2324$target_end_date[wday(df_hosp_2324$target_end_date) == '2'])
forecast_dates <- forecast_dates[as.Date('2023-01-01') < forecast_dates & 
                                   forecast_dates < as.Date('2024-05-01')]

naive_h2_quant_list <- list()  # Use a list to collect tibbles
i <- 1  # Counter for list indexing

for (this_loc in unique(df_hosp_versioned$geo_value)) {
  print(this_loc)
  for (this_date in forecast_dates) {
    this_date <- as.Date(this_date, origin = "1970-01-01")
    this_df <- df_hosp_versioned %>%
      filter(as_of == this_date, geo_value == this_loc,
             time_value < this_date, !is.na(value)) 
    if (nrow(this_df) < 3) next()
    nb_fit <- fit_nbinom(this_df$value, if_all_zeros='add1')
    this_quant <- tibble(
      outcome = 'hosp',
      location = name_map$location_number[which(name_map$geo_value == this_loc)],
      location_name = name_map$location_name[which(name_map$geo_value == this_loc)],
      forecast_date = this_date,
      horizon = these_target_days/7 - 1,
      target = paste0(these_target_days, ' day ahead inc hosp'),
      target_end_date = this_date + these_target_days,
      nbinom_mu = nb_fit$mu,
      nbinom_size = nb_fit$size
    ) %>%
      mutate(qv = map2(nbinom_mu, nbinom_size, ~ nbinom_to_quant(mu=.x, size=.y, quantiles=quantiles_h))) %>%
      unnest(qv) %>%
      dplyr::select(-nbinom_mu, -nbinom_size)
    
    naive_h2_quant_list[[i]] <- this_quant
    i <- i + 1
  }
}

naive_h2_quant <- bind_rows(naive_h2_quant_list)

#######################################
### hospitalizations 2024-2025
df_hosp_2425 <- read_parquet('surveillance data/versioned_data/covid_hosp_24-25.parquet') %>%
  mutate(obs_value = ifelse(value < 0, NA, value)) %>%
  filter(target_end_date <= '2024-05-04' | '2024-11-09' <= target_end_date) %>%
  arrange(target_end_date)

forecast_dates_2425 <- unique(df_hosp_2425$as_of) # Wednesday due dates
forecast_dates_2425 <- forecast_dates_2425[as.Date('2024-11-20') <= forecast_dates_2425]

quantiles_h <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
                 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

naive_2425_quant <- tibble()
for (this_loc in unique(df_hosp_2425$geo_value)) {
  print(this_loc)
  for (this_date in forecast_dates_2425) {
    this_date <- as.Date(this_date, origin = "1970-01-01")
    #this_loc = 'us'; this_date = forecast_dates_2425[1] 
    this_df <- filter(df_hosp_2425, geo_value == this_loc, as_of == this_date,
                      target_end_date < this_date, !is.na(value))
    if (nrow(this_df) < 3) next()
    nb_fit <- fit_nbinom(this_df$value, if_all_zeros='add1')
    this_quant <- tibble(
      outcome = 'hosp',
      location = name_map$location_number[which(name_map$geo_value == this_loc)],
      location_name = name_map$location_name[which(name_map$geo_value == this_loc)],
      forecast_date = this_date,
      reference_date = this_date + 3,
      horizon = (1:horizons_to_forecast) - 1,
      target = paste0((1:horizons_to_forecast) - 1, ' wk ahead inc covid hosp'),
      target_end_date = reference_date + ((1:horizons_to_forecast) - 1) * 7,
      nbinom_mu = nb_fit$mu,
      nbinom_size = nb_fit$size
    ) %>%
      mutate(qv = map2(nbinom_mu, nbinom_size, ~ nbinom_to_quant(mu=.x, size=.y, quantiles=quantiles_h))) %>%
      unnest(qv) %>%
      dplyr::select(-nbinom_mu, -nbinom_size)
    
    naive_2425_quant <- bind_rows(naive_2425_quant, this_quant)
  }
}


naive_covid_hosp_combined <- bind_rows(naive_h_quant, naive_h2_quant, naive_2425_quant)

if (horizons_to_forecast == 4) {
  write_parquet(naive_covid_hosp_combined, 'benchmark forecasts/naive/naive_covid_hosp.parquet')
} else {
  write_parquet(naive_covid_hosp_combined, 
    paste0('benchmark forecasts/naive/naive_covid_hosp_h', horizons_to_forecast, '.parquet'))
}


