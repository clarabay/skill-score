
library(arrow)
library(epidatr)
library(dplyr)
library(MMWRweek)
library(lubridate)

#### Get forecast dates (dates on which we want versioned data) ####
# Deaths
df_death <- read_parquet('surveillance data/observed data/truth_deaths.parquet') %>%
  mutate(obs_value = ifelse(value < 0, NA, value)) %>%
  arrange(target_end_date)

forecast_dates_covid_death <- unique(df_death$target_end_date) -5 

# Cases
df_case <- read_parquet('surveillance data/observed data/truth_cases.parquet') %>%
  mutate(obs_value = ifelse(value < 0, NA, value)) %>%
  arrange(target_end_date) %>%
  filter(geo_type == 'state')

forecast_dates_covid_cases <- unique(df_case$target_end_date) -5 

# Hospitalizations
df_hosp <- read_parquet('surveillance data/observed data/truth_hosp.parquet') %>%
  mutate(obs_value = ifelse(value < 0, NA, value)) %>%
  arrange(target_end_date)

forecast_dates_covid_hosp <- unique(df_hosp$target_end_date[wday(df_hosp$target_end_date) == '2'])

#df_hosp_2324 <- read_parquet('surveillance data/observed data/covid_hosp.parquet') %>%
  #rename(target_end_date = date) %>%
  #mutate(
  #  target_end_date = as.Date(target_end_date),
  #  obs_value = ifelse(value < 0, NA, value)) %>%
 # arrange(target_end_date)

df_hosp_2324 <- read_parquet('surveillance data/observed data/COVID19-observed.parquet') %>%
  filter(outcome == 'hosp') %>%
  mutate(
    target_end_date = as.Date(target_end_date),
    obs_value = ifelse(obs_value < 0, NA, value)) %>%
  arrange(target_end_date)

forecast_dates <- unique(df_hosp_2324$target_end_date[wday(df_hosp_2324$target_end_date) == '2'])
forecast_dates_covid_hosp_23_24 <- forecast_dates[as.Date('2023-01-01') < forecast_dates & 
                                   forecast_dates < as.Date('2024-05-01')]

rm(forecast_dates)

#### Request data from Delphi EpiData API ####
# Deaths
covid_sources <- covidcast_epidata()
covid_sources$signals$`jhu-csse:confirmed_incidence_num`

data_deaths <- pub_covidcast(
  source = 'jhu-csse',
  signals = "deaths_incidence_num",
  time_type = "day",
  # time_values = epirange("2020-08-08", "2021-02-28"),
  geo_type = "nation",
  as_of = "2020-04-02"
)

# Initialize an empty list
df_list <- list()

for (i in 1:length(forecast_dates_covid_death)) {
  print(i/length(forecast_dates_covid_death))
  date = as.Date(forecast_dates_covid_death[i])
  
  # Get nation data
  df_nation <- pub_covidcast(
    source = 'jhu-csse',
    signals = "deaths_incidence_num",
    time_type = "day",
    geo_type = "nation",
    as_of = date
  )
  df_nation$as_of <- date
  # Remove unnecessary columns
  df_nation <- df_nation[, !names(df_nation) %in% c("source", "geo_type", "direction", 
                                                    "lag", "stderr", "sample_size", "time_type")]
  
  # Get state data
  df_state <- pub_covidcast(
    source = 'jhu-csse',
    signals = "deaths_incidence_num",
    time_type = "day",
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
  
  # Add to list (much faster than rbinding to growing dataframe)
  df_list[[length(df_list) + 1]] <- this_df
}

# Combine all dataframes at once at the end
deaths <- do.call(rbind, df_list)

write_parquet(deaths, 'surveillance data/versioned_data/covid_deaths.parquet')
rm(deaths)


# Cases

data_cases <- pub_covidcast(
  source = 'jhu-csse',
  signals = "confirmed_incidence_num",
  time_type = "day",
  # time_values = epirange("2020-08-08", "2021-02-28"),
  geo_type = "nation",
  # as_of = "2020-04-02"
)


# Initialize an empty list
df_list <- list()

for (i in 1:length(forecast_dates_covid_cases)) {
  print(i/length(forecast_dates_covid_cases))
  date = as.Date(forecast_dates_covid_cases[i], origin = "1970-01-01")
  
  # Get nation data
  df_nation <- pub_covidcast(
    source = 'jhu-csse',
    signals = "confirmed_incidence_num",
    time_type = "day",
    geo_type = "nation",
    as_of = date
  )
  df_nation$as_of <- date
  # Remove unnecessary columns
  df_nation <- df_nation[, !names(df_nation) %in% c("source", "geo_type", "direction", 
                                                    "lag", "stderr", "sample_size", "time_type")]
  
  # Get state data
  df_state <- pub_covidcast(
    source = 'jhu-csse',
    signals = "confirmed_incidence_num",
    time_type = "day",
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
  
  # Add to list (much faster than rbinding to growing dataframe)
  df_list[[length(df_list) + 1]] <- this_df
}

# Combine all dataframes at once at the end
cases <- do.call(rbind, df_list)
cases$missing_sample_size <- cases$missing_stderr <- cases$missing_value <- NULL

write_parquet(cases, 'surveillance data/versioned_data/covid_cases.parquet')
rm(cases)

# Hospitalizations

covid_sources$signals$`nhsn:confirmed_admissions_covid_ew_prelim`
hosp <- pub_covidcast(
  source = 'hhs',
  signals = "confirmed_admissions_covid_1d",
  time_type = "day",
  # time_values = epirange("2020-10-01", "2025-04-15"),
  geo_type = "nation",
  as_of = "2020-11-16"
)

# Initialize an empty list
df_list <- list()

forecast_dates <- c(forecast_dates_covid_hosp, forecast_dates_covid_hosp_23_24)

for (i in 1:length(forecast_dates)) {
  print(i/length(forecast_dates))
  date = as.Date(forecast_dates[i], origin = "1970-01-01")
  
  # Get nation data
  df_nation <- pub_covidcast(
    source = 'hhs',
    signals = "confirmed_admissions_covid_1d",
    time_type = "day",
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
    signals = "confirmed_admissions_covid_1d",
    time_type = "day",
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
  
  # Add to list (much faster than rbinding to growing dataframe)
  df_list[[length(df_list) + 1]] <- this_df
}

# Combine all dataframes at once at the end
hosp <- do.call(rbind, df_list)
hosp$missing_sample_size <- hosp$missing_stderr <- hosp$missing_value <- NULL

write_parquet(hosp, 'surveillance data/versioned_data/covid_hosp.parquet')
rm(hosp)

# 2024-2025 Hospitalizations

# forecast_dates <- seq(as.Date("2024-11-30"), as.Date("2025-07-12"), by=7) - 3 # wednesdays
forecast_dates <- seq(as.Date("2024-10-12"), as.Date("2025-12-31"), by=7) - 3 # wednesdays
all_dates = list()
i = 1
for (d in forecast_dates) {
  date = as.Date(d, origin = "1970-01-01")
  date = paste0(epiyear(date), sprintf("%02d", epiweek(date)))
  nation <- pub_covidcast(
    source = 'nhsn',
    signals = "confirmed_admissions_covid_ew_prelim",
    time_type = "week",
    # time_values = epirange("2020-10-01", "2025-04-15"),
    geo_type = "nation",
    as_of = date
  )
  nation$as_of <- as.Date(d, origin = "1970-01-01")
  nation$geo_type <- nation$direction <- nation$lag <- nation$missing_sample_size <- nation$missing_stderr <- nation$missing_value <- nation$stderr <- nation$sample_size <- NULL
  state <- pub_covidcast(
    source = 'nhsn',
    signals = "confirmed_admissions_covid_ew_prelim",
    time_type = "week",
    # time_values = epirange("2020-10-01", "2025-04-15"),
    geo_type = "state",
    as_of = date
  )
  state$as_of <- as.Date(d, origin = "1970-01-01")
  state$geo_type <- state$direction <- state$lag <- state$missing_sample_size <- state$missing_stderr <- state$missing_value <- state$stderr <- state$sample_size <- NULL
  
  df = rbind(nation,state)
  all_dates[[i]] = df
  print(i/length(forecast_dates))
  i = i + 1
  # write.csv(df, paste0("versioned_data/covid_hosp25_", as.Date(d), '.csv'))
}
all_dates_df <- do.call(bind_rows, all_dates)
all_dates_df$target_end_date <- all_dates_df$time_value + 6
write_parquet(all_dates_df, "surveillance data/versioned_data/covid_hosp_24-25.parquet")

