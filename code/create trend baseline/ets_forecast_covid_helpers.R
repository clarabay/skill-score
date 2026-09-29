# Helpers for COVID ETS (additive damped): weekly deaths/cases, daily hospitalizations, weekly 24-25 hosp.
# Sourced from ets_additive_damped_forecast_covid.R

quantiles_covid <- c(
  0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
  0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99
)
forecast_levels <- 100 * c(1 - 2 * quantiles_covid[1:12])
STL_PERIOD_WEEKLY <- 52L
STL_PERIOD_DAILY <- 7L

stl_future_seasonal <- function(seas_vec, period, h_max) {
  last_cycle <- tail(seas_vec, period)
  idx <- ((seq_len(h_max) - 1L) %% period) + 1L
  as.numeric(last_cycle[idx])
}

stl_prefilter_log1p <- function(y, period, h_max) {
  n <- length(y)
  if (n < 2L * period) return(NULL)
  x <- log(y + 1)
  fit <- tryCatch(
    stats::stl(stats::ts(x, frequency = period), s.window = "periodic", robust = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  seas <- as.numeric(fit$time.series[, "seasonal"])
  list(
    x_deseasoned = x - seas,
    S_future = stl_future_seasonal(seas, period, h_max)
  )
}

stl_prefilter_level_add <- function(y, period, h_max) {
  n <- length(y)
  if (n < 2L * period) return(NULL)
  y_pos <- y + 1
  fit <- tryCatch(
    stats::stl(stats::ts(y_pos, frequency = period), s.window = "periodic", robust = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  seas <- as.numeric(fit$time.series[, "seasonal"])
  list(
    x_deseasoned = as.numeric(y_pos - seas),
    S_future = stl_future_seasonal(seas, period, h_max)
  )
}

# Weekly targets: 1–h_max wk ahead; same week-ending convention as naive_covid (base = forecast_date + 5).
pred_to_quant_covid_weekly <- function(pred, meta, log_scale, forecast_adjustment, h_max, inc_suffix) {
  if (!is.null(forecast_adjustment)) {
    if (length(forecast_adjustment) != h_max) {
      stop("forecast_adjustment length must match h_max", call. = FALSE)
    }
    pred$mean <- as.numeric(pred$mean) + forecast_adjustment
    if (!is.null(pred$lower)) pred$lower <- sweep(as.matrix(pred$lower), 1L, forecast_adjustment, "+")
    if (!is.null(pred$upper)) pred$upper <- sweep(as.matrix(pred$upper), 1L, forecast_adjustment, "+")
  }
  out <- bind_rows(
    as_tibble(pred$lower) %>%
      mutate(horizon = 0:(h_max - 1), lower = TRUE) %>%
      pivot_longer(cols = !c(horizon, lower)),
    as_tibble(pred$upper) %>%
      mutate(horizon = 0:(h_max - 1), lower = FALSE) %>%
      pivot_longer(cols = !c(horizon, lower))
  ) %>%
    rename(quantile = name) %>%
    filter(!(quantile == "0%" & lower)) %>%
    mutate(
      location = meta$location,
      location_name = meta$location_name,
      forecast_date = meta$forecast_date,
      reference_date = meta$reference_date,
      outcome = meta$outcome,
      target = paste0(horizon + 1, " wk ahead inc ", inc_suffix),
      target_end_date = meta$target_end_date_base + horizon * 7,
      quantile = as.numeric(str_remove(quantile, "%")) / 100,
      quantile = ifelse(lower, (1 - quantile) / 2, 1 - (1 - quantile) / 2)
    )
  if (log_scale) {
    out <- out %>% mutate(value = ifelse(lower, floor(exp(value) - 1), ceiling(exp(value) - 1)))
  } else {
    out <- out %>% mutate(value = ifelse(lower, floor(value), ceiling(value)))
  }
  out %>%
    mutate(value = pmax(0, value)) %>%
    arrange(horizon, quantile) %>%
    dplyr::select(-lower, -horizon)
}

# Daily hosp: horizons 0..h_max-1 map to target_end_date = forecast_date + horizon + 1.
pred_to_quant_covid_daily_hosp <- function(pred, meta, log_scale, forecast_adjustment, h_max) {
  if (!is.null(forecast_adjustment)) {
    if (length(forecast_adjustment) != h_max) {
      stop("forecast_adjustment length must match h_max", call. = FALSE)
    }
    pred$mean <- as.numeric(pred$mean) + forecast_adjustment
    if (!is.null(pred$lower)) pred$lower <- sweep(as.matrix(pred$lower), 1L, forecast_adjustment, "+")
    if (!is.null(pred$upper)) pred$upper <- sweep(as.matrix(pred$upper), 1L, forecast_adjustment, "+")
  }
  fd <- meta$forecast_date
  out <- bind_rows(
    as_tibble(pred$lower) %>%
      mutate(horizon = 0:(h_max - 1), lower = TRUE) %>%
      pivot_longer(cols = !c(horizon, lower)),
    as_tibble(pred$upper) %>%
      mutate(horizon = 0:(h_max - 1), lower = FALSE) %>%
      pivot_longer(cols = !c(horizon, lower))
  ) %>%
    rename(quantile = name) %>%
    filter(!(quantile == "0%" & lower)) %>%
    mutate(
      location = meta$location,
      location_name = meta$location_name,
      forecast_date = fd,
      reference_date = meta$reference_date,
      outcome = meta$outcome,
      target = paste0(horizon + 1, " day ahead inc hosp"),
      target_end_date = as.Date(fd) + horizon + 1L,
      quantile = as.numeric(str_remove(quantile, "%")) / 100,
      quantile = ifelse(lower, (1 - quantile) / 2, 1 - (1 - quantile) / 2)
    )
  if (log_scale) {
    out <- out %>% mutate(value = ifelse(lower, floor(exp(value) - 1), ceiling(exp(value) - 1)))
  } else {
    out <- out %>% mutate(value = ifelse(lower, floor(value), ceiling(value)))
  }
  out %>%
    mutate(value = pmax(0, value)) %>%
    arrange(horizon, quantile) %>%
    dplyr::select(-lower, -horizon)
}

# 2024–25 weekly COVID hosp: targets "0"–"3" wk ahead inc covid hosp; reference_date = forecast_date + 3.
pred_to_quant_covid_hosp2425 <- function(pred, meta, log_scale, forecast_adjustment, h_max) {
  if (!is.null(forecast_adjustment)) {
    if (length(forecast_adjustment) != h_max) {
      stop("forecast_adjustment length must match h_max", call. = FALSE)
    }
    pred$mean <- as.numeric(pred$mean) + forecast_adjustment
    if (!is.null(pred$lower)) pred$lower <- sweep(as.matrix(pred$lower), 1L, forecast_adjustment, "+")
    if (!is.null(pred$upper)) pred$upper <- sweep(as.matrix(pred$upper), 1L, forecast_adjustment, "+")
  }
  ref <- meta$reference_date
  out <- bind_rows(
    as_tibble(pred$lower) %>%
      mutate(horizon = 0:(h_max - 1), lower = TRUE) %>%
      pivot_longer(cols = !c(horizon, lower)),
    as_tibble(pred$upper) %>%
      mutate(horizon = 0:(h_max - 1), lower = FALSE) %>%
      pivot_longer(cols = !c(horizon, lower))
  ) %>%
    rename(quantile = name) %>%
    filter(!(quantile == "0%" & lower)) %>%
    mutate(
      location = meta$location,
      location_name = meta$location_name,
      forecast_date = meta$forecast_date,
      reference_date = ref,
      outcome = meta$outcome,
      target = paste0(horizon, " wk ahead inc covid hosp"),
      target_end_date = as.Date(ref) + horizon * 7,
      quantile = as.numeric(str_remove(quantile, "%")) / 100,
      quantile = ifelse(lower, (1 - quantile) / 2, 1 - (1 - quantile) / 2)
    )
  if (log_scale) {
    out <- out %>% mutate(value = ifelse(lower, floor(exp(value) - 1), ceiling(exp(value) - 1)))
  } else {
    out <- out %>% mutate(value = ifelse(lower, floor(value), ceiling(value)))
  }
  out %>%
    mutate(value = pmax(0, value)) %>%
    arrange(horizon, quantile) %>%
    dplyr::select(-lower, -horizon)
}

forecast_additive_covid_one <- function(y_raw, meta, use_stl, use_log, h_max, stl_period, pred_fn) {
  if (length(y_raw) < 4L) return(NULL)
  # Daily hosp: STL/ETS cost grows sharply with series length; fit on recent history only.
  if (identical(stl_period, STL_PERIOD_DAILY) && length(y_raw) > 800L) {
    y_raw <- tail(y_raw, 800L)
  }
  stl_ok <- isTRUE(use_stl) && length(y_raw) >= 2L * stl_period

  if (use_log) {
    stl_log <- if (stl_ok) stl_prefilter_log1p(y_raw, stl_period, h_max) else NULL
    x <- log(y_raw + 1)
    if (!is.null(stl_log)) x <- stl_log$x_deseasoned
    adj <- if (!is.null(stl_log)) stl_log$S_future else rep(0, h_max)
    if (all(x == 0)) x <- c(x, log(1 + 1))
    if (any(!is.finite(x))) return(NULL)
    m <- tryCatch(ets(x, model = "AAN", damped = TRUE), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    pred <- tryCatch(
      forecast(m, h = h_max, level = forecast_levels),
      error = function(e) NULL
    )
    if (is.null(pred)) return(NULL)
    pred_fn(pred, meta, TRUE, adj)
  } else {
    stl_lev <- if (stl_ok) stl_prefilter_level_add(y_raw, stl_period, h_max) else NULL
    x <- y_raw + 1
    if (!is.null(stl_lev)) x <- stl_lev$x_deseasoned
    adj <- if (!is.null(stl_lev)) stl_lev$S_future else rep(0, h_max)
    if (any(!is.finite(x))) return(NULL)
    m <- tryCatch(ets(x, model = "AAN", damped = TRUE), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    pred <- tryCatch(
      forecast(m, h = h_max, level = forecast_levels),
      error = function(e) NULL
    )
    if (is.null(pred)) return(NULL)
    pred_fn(pred, meta, FALSE, adj)
  }
}

# Log-scale ETS point forecast only (no PI). Uses full `y_raw` for all streams — no daily tail(800)
# truncation (see `forecast_additive_covid_one` for the performance-oriented variant).
# Used by `ets_additive_damped_empirical_log_covid.R` for empirical PI construction.
ets_covid_log_point_mean <- function(y_raw, h_max) {
  if (length(y_raw) < 4L) return(NULL)
  x <- log(y_raw + 1)
  if (all(x == 0)) x <- c(x, log(1 + 1))
  if (any(!is.finite(x))) return(NULL)
  m <- tryCatch(ets(x, model = "AAN", damped = TRUE), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  pred <- tryCatch(forecast(m, h = h_max), error = function(e) NULL)
  if (is.null(pred)) return(NULL)
  as.numeric(pred$mean)
}
