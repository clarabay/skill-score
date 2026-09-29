rm(list = ls())

require(tidyverse)
require(arrow)

source("code/WIS.R")

# ── 1. Observed data (COVID-19 cases, state-level only) ─────────────────────
df_obs <- read_parquet("surveillance data/observed data/COVID19-observed.parquet") %>%
  filter(outcome == "case", nchar(location) == 2L) %>%
  arrange(target_end_date)

# ── 2. Hindcast WIS (keyed by location × outcome × target_end_date) ──────────
hindcast <- read_parquet("benchmark forecasts/hindcasts/COVID19-case-state-hindcast-quantiles.parquet") %>%
  mutate(log_value = log(value + 1), outcome = "case") %>%
  arrange(target_end_date)

wis_hindcast <- hindcast %>%
  left_join(
    df_obs %>% dplyr::select(location, target_end_date, obs_value, log_obs_value),
    by = c("location", "target_end_date")
  ) %>%
  group_by(location, outcome, target_end_date) %>%
  summarize(
    wis_hindcast = weighted_interval_score(quantile, value, obs_value),
    log_wis_hindcast = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups = "drop"
  )

# ── 3. Naive WIS (from combined parquet for consistency across all models) ─────
wis_naive <- read_parquet("scored forecasts/COVID19_subset_combined_WIS.parquet") %>%
  filter(outcome == "case", !is.na(log_wis_naive)) %>%
  distinct(location, outcome, forecast_date, target_end_date, .keep_all = TRUE) %>%
  dplyr::select(location, outcome, forecast_date, target_end_date, wis_naive, log_wis_naive)

# ── 4. Log-baseline WIS (COVID cases) ───────────────────────────────────────
log_baseline <- read_parquet("benchmark forecasts/random walk baseline/covid_case_log_baseline_h12.pq") %>%
  mutate(
    log_value = log(value + 1),
    forecast_date = as.Date(forecast_date),
    target_end_date = as.Date(target_end_date),
    outcome = "case"
  ) %>%
  filter(horizon %in% 1:4) %>%
  arrange(target_end_date)

wis_log_baseline <- log_baseline %>%
  left_join(
    df_obs %>% dplyr::select(location, target_end_date, obs_value, log_obs_value),
    by = c("location", "target_end_date")
  ) %>%
  group_by(location, outcome, forecast_date, target_end_date) %>%
  summarize(
    wis_lbaseline = weighted_interval_score(quantile, value, obs_value),
    log_wis_lbaseline = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups = "drop"
  )

# ── 5a. Score Log-baseline model (as forecast to plot) ──────────────────────
log_baseline_horizon <- log_baseline %>%
  mutate(horizon = as.numeric(horizon) - 1L)

wis_log_baseline_model <- log_baseline_horizon %>%
  left_join(
    df_obs %>% dplyr::select(location, target_end_date, obs_value, log_obs_value),
    by = c("location", "target_end_date")
  ) %>%
  group_by(location, outcome, forecast_date, target_end_date, horizon) %>%
  summarize(
    location_name = first(location_name),
    target = first(target),
    wis = weighted_interval_score(quantile, value, obs_value),
    log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
    .groups = "drop"
  ) %>%
  mutate(
    team = "Log-baseline",
    reference_date = as.Date(forecast_date)
  )

# ── 5b. Score the 19 cBase models ──────────────────────────────────────────
cbase_dir <- "baseline_models"
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

cbase_wis_list <- list()

cat("Scoring cBase models...\n")
for (mid in model_ids) {
  fname <- file.path(cbase_dir, paste0("cbase_model_", mid, ".csv"))
  cat("  model", mid, ":", model_names[as.character(mid)], "\n")

  raw <- read_csv(fname, show_col_types = FALSE) %>%
    mutate(
      forecast_date = as.Date(forecast_date),
      target_end_date = as.Date(target_end_date),
      reference_date = as.Date(reference_date)
    )

  scored <- raw %>%
    left_join(
      df_obs %>% dplyr::select(location, target_end_date, obs_value, log_obs_value),
      by = c("location", "target_end_date")
    ) %>%
    group_by(location, location_name, reference_date, forecast_date, target_end_date, horizon) %>%
    summarize(
      target = first(target),
      wis = weighted_interval_score(quantile, value, obs_value),
      log_wis = weighted_interval_score(quantile, log_value, log_obs_value),
      .groups = "drop"
    ) %>%
    mutate(
      team = model_names[as.character(mid)],
      outcome = "case"
    )

  cbase_wis_list[[mid]] <- scored
}
cat("Done scoring.\n")

wis_cbase <- bind_rows(cbase_wis_list)

# ── 6. Merge reference scores onto all forecast models ──────────────────────
cbase_combined <- wis_cbase %>%
  left_join(wis_hindcast, by = c("location", "outcome", "target_end_date")) %>%
  left_join(wis_naive, by = c("location", "outcome", "forecast_date", "target_end_date")) %>%
  left_join(wis_log_baseline, by = c("location", "outcome", "forecast_date", "target_end_date"))

log_baseline_combined <- wis_log_baseline_model %>%
  left_join(wis_hindcast, by = c("location", "outcome", "target_end_date")) %>%
  left_join(wis_naive, by = c("location", "outcome", "forecast_date", "target_end_date")) %>%
  left_join(wis_log_baseline, by = c("location", "outcome", "forecast_date", "target_end_date"))

# ── 7. Load existing combined scores (COVIDhub ensemble + baseline) ────────
existing_combined <- read_parquet("scored forecasts/COVID19_subset_combined_WIS.parquet") %>%
  filter(outcome == "case", team %in% c("COVIDhub-4_week_ensemble", "COVIDhub-baseline"))

# ── 8. Stack everything together ────────────────────────────────────────────
shared_cols <- c(
  "team", "location", "location_name", "outcome", "forecast_date", "reference_date",
  "target_end_date", "target", "horizon",
  "wis", "log_wis",
  "wis_hindcast", "log_wis_hindcast",
  "wis_naive", "log_wis_naive",
  "wis_lbaseline", "log_wis_lbaseline"
)

existing_slim <- existing_combined %>%
  dplyr::select(any_of(shared_cols))

cbase_slim <- cbase_combined %>%
  dplyr::select(any_of(shared_cols))

log_baseline_slim <- log_baseline_combined %>%
  dplyr::select(any_of(shared_cols))

all_models <- bind_rows(existing_slim, log_baseline_slim, cbase_slim)

# ── 9. Assign seasons & filter (COVID case challenge dates) ─────────────────
all_models <- all_models %>%
  mutate(season = format(target_end_date, "%Y")) %>%
  filter(
    !is.na(forecast_date),
    "2020-07-28" <= forecast_date,
    forecast_date <= "2021-12-21"
  ) %>%
  filter(season %in% c("2020", "2021")) %>%
  filter(horizon %in% 0:3) %>%
  group_by(team, location, target_end_date) %>%
  mutate(n_horizon = n()) %>%
  ungroup() %>%
  filter(n_horizon == 4L) %>%
  dplyr::select(-n_horizon)

# ── 10. Save scored results ─────────────────────────────────────────────────
write_parquet(all_models, "scored forecasts/cBase_combined_WIS.parquet")
cat("Saved: scored forecasts/cBase_combined_WIS.parquet\n")

cat("\nScoring complete. Run 'code/Plot cBase models.R' for plots.\n")
