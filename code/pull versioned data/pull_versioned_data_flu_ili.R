####

library(arrow)
library(epidatr)
library(dplyr)
library(MMWRweek)
library(lubridate)


# Import
source('code/delphi_epidata.R') 
# Fetch data

forecast_dates <- bind_rows(
  tibble(
    year = 2015,
    epi_week = 42:52),
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


forecast_dates$epiyearweek <- as.integer(paste0(forecast_dates$year, sprintf("%02d", forecast_dates$epi_week)))

regions = list('nat',
               'hhs1',
               'hhs2',
               'hhs3',
               'hhs4',
               'hhs5',
               'hhs6',
               'hhs7',
               'hhs8',
               'hhs9',
               'hhs10'
)


res <- Epidata$fluview(regions, list(Epidata$range(199740, 202013)),
                       issues = list(forecast_dates$epiyearweek), auth = '2baab43442a82')

# res <- Epidata$fluview(regions, list(201445),
#                        issues = list(201543), auth = '2baab43442a82')
# res$epidata[[1]]


cat(paste(res$result, res$message, length(res$epidata), "\n"))

df_versioned <- bind_rows(res$epidata)
df_versioned$lag <- NULL
df_versioned <- arrange(df_versioned, issue)

df_versioned$num_age_0 <- df_versioned$num_age_1 <- df_versioned$num_age_2 <- df_versioned$num_age_3 <- df_versioned$num_age_4 <- df_versioned$num_age_5 <- NULL

res <- Epidata$fluview(regions, list(Epidata$range(199740, 202013)),
                       auth = '2baab43442a82')
res$epidata[[1]]
df_latest <- bind_rows(res$epidata)
df_latest$lag <- NULL
df_latest$num_age_0 <- df_latest$num_age_1 <- df_latest$num_age_2 <- df_latest$num_age_4 <- df_latest$num_age_5 <- NULL

unique_issues <- sort(unique(df_versioned$issue))

all_versioned_list <- list()
for (this_ui in unique_issues) {
  issue <- as.numeric(issue)
  this_issue <- filter(df_versioned, issue==this_ui)
  min_date <- min(this_issue$epiweek)
  remainder <- filter(df_latest, epiweek < min_date)
  remainder$as_of <- this_ui
  this_issue$as_of <- this_ui
  all_versioned_list[[i]] <- bind_rows(this_issue, remainder)
  i = i+1
}

all_versioned_df <- do.call(bind_rows, all_versioned_list) %>% arrange(as_of, epiweek)

write_parquet(all_versioned_df, "surveillance data/versioned_data/influenza_ili_new.parquet")
