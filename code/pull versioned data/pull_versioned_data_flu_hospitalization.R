library(arrow)
library(epidatr)
library(dplyr)
library(MMWRweek)
library(lubridate)

# INFLUENZA
df_obs <- bind_rows( # combine rows from both parquet files
  read_csv('surveillance data/observed data/FluSight-hosp-2022-2023-observed.csv') %>% # read 
    filter(target_end_date < as.Date('2023-06-01')) %>% # ignore dates after 2023-06-01
    mutate(forecast_date = target_end_date - 5), # forecast date is 5 days before end of epiweek, so Monday
  read_csv('surveillance data/observed data/FluSight-hosp-2025-2026-observed.csv') %>% # read
    mutate(
      target_end_date = as.Date(target_end_date), 
      forecast_date = target_end_date - 3) %>% # forecast date is 3 days before end of epiweek, so Wednesday
    filter(as.Date('2023-06-01') <= target_end_date) # ignore dates before 2023-06-01
) %>% # take combined dataframe
  mutate(
    obs_value = ifelse(value < 0, NA, value), # convert negative values to NA
    log_obs_value = log(obs_value + 1)) %>% # create log transformed column
  dplyr::select(-value) %>%
  arrange(target_end_date) # target end dates are saturdays 

forecast_dates <- unique(df_obs$forecast_date)
forecast_dates <- forecast_dates[
  (as.Date('2022-01-01') <= forecast_dates & forecast_dates < as.Date('2023-07-01')) |
    (as.Date('2023-10-01') <= forecast_dates & forecast_dates < as.Date('2024-06-01'))]
forecast_dates

# data <- pub_covidcast(
#   source = 'hhs',
#   signals = "confirmed_admissions_influenza_1d",
#   time_type = "day",
#   # time_values = epirange("2020-10-01", "2021-09-25"),
#   geo_type = "nation",
#   as_of = "2020-11-16"
# )

# Initialize an empty list
df_list <- list()

# forecast_dates <- c(forecast_dates_covid_hosp, forecast_dates_covid_hosp_23_24)

# from_date = as.Date("2020-08-08")
# from_date = paste0(epiyear(from_date), sprintf("%02d", epiweek(from_date)))
# to_date = as.Date("2022-06-06")
# to_date = paste0(epiyear(to_date), sprintf("%02d", epiweek(to_date)))
# 
from_date = "2020-08-08"
date = forecast_dates[1]
df_ca <- pub_covidcast(
  source = 'hhs',
  signals = "confirmed_admissions_influenza_1d",
  time_type = "day",
  time_values = epirange(from_date, date),
  geo_type = "state",
  geo_values = 'ca',
  as_of = date
)
df_ca <- df_ca %>%
  mutate(value = ifelse(value < 0, NA, value)) %>%
  mutate(target_end_date = ceiling_date(time_value, "week", week_start = 7)) %>%
  mutate(target_end_date = target_end_date-1) %>%
  group_by(geo_value, target_end_date) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = 'drop') %>%
  arrange(geo_value)

time_type = 'day'
for (i in 1:length(forecast_dates)) {
  print(i/length(forecast_dates))
  date = as.Date(forecast_dates[i], origin = "1970-01-01")
  from_date = "2020-08-08"
  # Get nation data
  df_nation <- pub_covidcast(
    source = 'hhs',
    signals = "confirmed_admissions_influenza_1d",
    time_type = time_type,
    time_values = epirange(from_date, date),
    geo_type = "nation",
    as_of = date
  )
  df_nation$as_of <- date
  # Remove unnecessary columns
  df_nation <- df_nation[, !names(df_nation) %in% c("source", "geo_type", "direction", 
                                                    "lag", "stderr", "sample_size", "time_type")]
  
  # Get state data
  df_state <- pub_covidcast(
    source = 'hhs',
    signals = "confirmed_admissions_influenza_1d",
    time_type = time_type,
    time_values = epirange(from_date, date),
    geo_type = "state",
    as_of = date
  )
  df_state$as_of <- date
  # Remove unnecessary columns
  df_state <- df_state[, !names(df_state) %in% c("source", "geo_type", "direction", 
                                                 "lag", "stderr", "sample_size", "time_type")]
  
  # Combine state and nation for this date
  this_df <- rbind(df_state, df_nation)
  this_df$as_of <- date
  
  # Add to list 
  df_list[[length(df_list) + 1]] <- this_df
}

# Combine all dataframes 
hosp <- do.call(rbind, df_list)
hosp$missing_sample_size <- hosp$missing_stderr <- hosp$missing_value <- NULL

# hosp = read_parquet("surveillance data/versioned_data/influenza_hosp.parquet")
# hosp <- hosp %>%
#   mutate(value = ifelse(value < 0, NA, value)) %>%
#   mutate(target_end_date = ceiling_date(time_value, "week", week_start = 7)) %>%
#   mutate(target_end_date = target_end_date-1) %>%
#   group_by(as_of, geo_value, target_end_date) %>%
#   summarise(value = sum(value, na.rm = TRUE), .groups = 'drop') %>%
#   arrange(geo_value)

write_parquet(hosp, 'surveillance data/versioned_data/influenza_hosp.parquet')
rm(hosp)

# 2024-2026
epiweek("2024-11-23")
# Reference dates (Saturdays) span both seasons: 2024-25 runs 2024-11-23 to 2025-05-31,
# 2025-26 runs 2025-11-22 to 2026-05-30 (per dat/flu24-25_ and flu25-26_number-models-submitted-weekly.csv,
# cross-checked against the flu_forecasts_24-25 / flu_forecasts_25-26 submission files).
# Minus 3 gives the Wednesday issue dates. Yields 80 as_of values, 2024-11-20 .. 2026-05-27.
forecast_dates <- seq(as.Date("2024-11-23"), as.Date("2026-05-30"), by=7) - 3 # wednesdays
all_dates = list()
i = 1
for (d in forecast_dates) {
  date = as.Date(d, origin = "1970-01-01")
  date = paste0(epiyear(date), sprintf("%02d", epiweek(date)))
  nation <- pub_covidcast(
    source = 'nhsn',
    signals = "confirmed_admissions_flu_ew_prelim",
    time_type = "week",
    # time_values = epirange("2020-10-01", "2025-04-15"),
    geo_type = "nation",
    as_of = date
  )
  nation$as_of <- as.Date(d, origin = "1970-01-01")
  nation$geo_type <- nation$direction <- nation$lag <- nation$missing_sample_size <- nation$missing_stderr <- nation$missing_value <- nation$stderr <- nation$sample_size <- NULL
  state <- pub_covidcast(
    source = 'nhsn',
    signals = "confirmed_admissions_flu_ew_prelim",
    time_type = "week",
    # time_values = epirange("2020-10-01", "2025-04-15"),
    geo_type = "state",
    as_of = date
  )
  state$as_of <- as.Date(d, origin = "1970-01-01")
  state$geo_type <- state$direction <- state$lag <- state$missing_sample_size <- state$missing_stderr <- state$missing_value <- state$stderr <- state$sample_size <- NULL
  
  df = rbind(nation,state)
  all_dates[[i]] = df
  i = i + 1
  # write.csv(df, paste0("surveillance data/versioned_data/covid_hosp25_", as.Date(d), '.csv'))
}
all_dates_df <- do.call(bind_rows, all_dates)
write_parquet(all_dates_df, "surveillance data/versioned_data/influenza_hosp_24-26.parquet")
