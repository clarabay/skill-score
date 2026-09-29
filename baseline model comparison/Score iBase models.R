rm(list = ls())

require(tidyverse)
require(arrow)

source("code/WIS.R")

# ── 1. Observed data ────────────────────────────────────────────────────────
df_obs <- bind_rows(
  read_csv("surveillance data/observed data/FluSight-hosp-2022-2023-observed.csv"),
  read_csv("surveillance data/observed data/FluSight-hosp-observed.csv")
) %>%
  mutate(
    obs_value = ifelse(value < 0, NA, value),
    log_obs_value = log(obs_value + 1)
  ) %>%
  dplyr::select(-value) %>%
  arrange(target_end_date)

# ── 2. Hindcast WIS (keyed by location × target_end_date) ──────────────────
hindcast <- read_parquet("benchmark forecasts/hindcasts/FluSight-hosp-hindcast-quantiles.parquet") %>%
  mutate(
    location_name = ifelse(location_name == "United States", "US", location_name),
    log_value = log(value + 1)
  ) %>%
  arrange(target_end_date)

wis_hindcast <- hindcast %>%
  left_join(df_obs, by = c("location", "target_end_date")) %>%
  group_by(location, target_end_date) %>%
  summarize(
    wis_hindcast = weighted_interval_score(quantile, value, obs_value),
    log_wis_hindcast = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups = "drop"
  )

# ── 3. Naive WIS (from combined parquet for consistency across all models) ────
wis_naive <- read_parquet("scored forecasts/FluSight-hosp_combined_WIS.parquet") %>%
  filter(!is.na(log_wis_naive)) %>%
  distinct(location, reference_date, target_end_date, .keep_all = TRUE) %>%
  dplyr::select(location, reference_date, target_end_date, wis_naive, log_wis_naive)

# ── 4. Log-baseline WIS ────────────────────────────────────────────────────
log_baseline <- read_parquet("benchmark forecasts/random walk baseline/flu_hosp_log_baseline_h12.pq") %>%
  mutate(
    log_value = log(value + 1),
    forecast_date = as.Date(forecast_date),
    reference_date = as.Date(reference_date),
    reference_date = ifelse(
      is.na(reference_date),
      ceiling_date(forecast_date, unit = "week", week_start = 6),
      reference_date
    ),
    reference_date = as.Date(reference_date)
  ) %>%
  arrange(target_end_date)

wis_log_baseline <- log_baseline %>%
  left_join(df_obs, by = c("location", "target_end_date")) %>%
  group_by(location, reference_date, target_end_date) %>%
  summarize(
    wis_lbaseline = weighted_interval_score(quantile, value, obs_value),
    log_wis_lbaseline = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups = "drop"
  )

# ── 5a. Score Log-baseline model ───────────────────────────────────────────
# Horizon is 1-12; use 1-4 and map to 0-3 for consistency with other models.
log_baseline_horizon <- log_baseline %>%
  filter(horizon %in% 1:4) %>%
  mutate(horizon = as.numeric(horizon) - 1L)

wis_log_baseline_model <- log_baseline_horizon %>%
  left_join(df_obs %>% dplyr::select(location, target_end_date, obs_value, log_obs_value),
            by = c("location", "target_end_date")) %>%
  group_by(location, reference_date, forecast_date, target_end_date, horizon) %>%
  summarize(
    location_name = first(location_name),
    target = first(target),
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups = "drop"
  ) %>%
  mutate(team = "Log-baseline")

# ── 5b. Score the 19 iBase models ───────────────────────────────────────────
ibase_dir <- "baseline_models"
model_ids <- 1:19

model_names <- c(
  "1"  = "Constant",
  "2"  = "Marginal",
  "3"  = "KDE",
  "4"  = "LSD(52,5)",
  "5"  = "OLS(5,1)",
  "6"  = "IDS(3)",
  "7"  = "ARMA(1,1)",
  "8"  = "INARCH(1)",
  "9"  = "ETS+ST",
  "10" = "STL+ST",
  "11" = "log-Constant",
  "12" = "log-Marginal",
  "13" = "log-KDE",
  "14" = "log-LSD(52,5)",
  "15" = "log-OLS(5,1)",
  "16" = "log-IDS(3)",
  "17" = "log-ARMA(1,1)",
  "18" = "log-ETS+ST",
  "19" = "log-STL+ST"
)

ibase_wis_list <- list()

cat("Scoring iBase models...\n")
for (mid in model_ids) {
  fname <- file.path(ibase_dir, paste0("ibase_model_", mid, ".csv"))
  cat("  model", mid, ":", model_names[as.character(mid)], "\n")

  raw <- read_csv(fname, show_col_types = FALSE) %>%
    mutate(
      forecast_date = as.Date(forecast_date),
      target_end_date = as.Date(target_end_date),
      reference_date = as.Date(reference_date)
    )

  scored <- raw %>%
    left_join(df_obs, by = c("location", "target_end_date")) %>%
    group_by(location, reference_date, forecast_date, target_end_date, horizon) %>%
    summarize(
      location_name = first(location_name.x),
      target = first(target),
      wis = weighted_interval_score(quantile, value, obs_value),
      log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
      .groups = "drop"
    ) %>%
    mutate(team = model_names[as.character(mid)])

  ibase_wis_list[[mid]] <- scored
}
cat("Done scoring.\n")

wis_ibase <- bind_rows(ibase_wis_list)

# ── 6. Merge reference scores onto all forecast models ───────────────────────
ibase_combined <- wis_ibase %>%
  left_join(wis_hindcast, by = c("location", "target_end_date")) %>%
  left_join(wis_naive, by = c("location", "reference_date", "target_end_date")) %>%
  left_join(wis_log_baseline, by = c("location", "reference_date", "target_end_date"))

log_baseline_combined <- wis_log_baseline_model %>%
  left_join(wis_hindcast, by = c("location", "target_end_date")) %>%
  left_join(wis_naive, by = c("location", "reference_date", "target_end_date")) %>%
  left_join(wis_log_baseline, by = c("location", "reference_date", "target_end_date"))

# ── 7. Load existing combined scores (ensemble + official baseline) ────────
existing_combined <- read_parquet("scored forecasts/FluSight-hosp_combined_WIS.parquet") %>%
  filter(team %in% c("FluSight-ensemble", "FluSight-baseline"))

# ── 8. Stack everything together ───────────────────────────────────────────
shared_cols <- c(
  "team", "location", "location_name", "forecast_date", "reference_date",
  "target_end_date", "target", "horizon",
  "wis", "log_wis",
  "wis_hindcast", "log_wis_hindcast",
  "wis_naive", "log_wis_naive",
  "wis_lbaseline", "log_wis_lbaseline"
)

existing_slim <- existing_combined %>%
  dplyr::select(any_of(shared_cols))

ibase_slim <- ibase_combined %>%
  dplyr::select(any_of(shared_cols))

log_baseline_slim <- log_baseline_combined %>%
  dplyr::select(any_of(shared_cols))

all_models <- bind_rows(existing_slim, log_baseline_slim, ibase_slim)

# ── 9. Assign seasons & filter (same logic as Score FluSight Hospitalizations.R)
all_models <- all_models %>%
  mutate(
    season = ifelse("2022-01-10" <= forecast_date & forecast_date <= "2022-06-20",
      "2021-2022", NA
    ),
    season = ifelse("2022-10-17" <= forecast_date & forecast_date <= "2023-05-17",
      "2022-2023", season
    ),
    season = ifelse("2023-10-01" <= reference_date & target_end_date < "2024-05-04",
      "2023-2024", season
    ),
    season = ifelse("2024-11-20" <= reference_date & reference_date <= "2025-05-31",
      "2024-2025", season
    )
  ) %>%
  filter(!is.na(season)) %>%
  filter(horizon %in% 0:3) %>%
  group_by(team, location, target_end_date) %>%
  mutate(n_horizon = n()) %>%
  ungroup() %>%
  filter(n_horizon == 4) %>%
  dplyr::select(-n_horizon) %>%
  filter(!(reference_date %in% as.Date(c("2022-01-15", "2022-01-22", "2025-05-31", "2023-10-07", "2025-01-25"))))

# ── 10. Save scored results ───────────────────────────────────────────────
write_parquet(all_models, "scored forecasts/iBase_combined_WIS.parquet")
cat("Saved: scored forecasts/iBase_combined_WIS.parquet\n")

cat("\nScoring complete. Run 'code/Plot iBase models.R' for plots.\n")
