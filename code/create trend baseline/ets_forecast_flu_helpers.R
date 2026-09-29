# Shared helpers for FluSight hospitalization ETS pipelines (additive / multiplicative).
# Sourced from ets_additive_damped_forecast_flu.R and ets_multiplicative_damped_forecast_flu.R.
# Run those scripts from the repository root.

horizons_to_forecast <- 12L
quantiles_flu <- c(
  0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
  0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99
)
forecast_levels <- 100 * c(1 - 2 * quantiles_flu[1:12])
STL_PERIOD <- 52L

stl_future_seasonal <- function(seas_vec, period, h_max) {
  last_cycle <- tail(seas_vec, period)
  idx <- ((seq_len(h_max) - 1L) %% period) + 1L
  as.numeric(last_cycle[idx])
}

# STL on log(y+1); additive deseason in log space; reattach seasonal in log space before expm1.
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

# STL on (y+1) in levels; additive seasonal for recomposition (pred + S_future).
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

# forecast::ets PI -> FluSight quantiles (matches naiveETS_forecast_flu.R).
# log_scale: TRUE if pred mean/intervals are on log(y+1) scale.
pred_to_quant <- function(pred, meta, log_scale = TRUE, forecast_adjustment = NULL) {
  if (!is.null(forecast_adjustment)) {
    if (length(forecast_adjustment) != horizons_to_forecast) {
      stop("forecast_adjustment must have length horizons_to_forecast", call. = FALSE)
    }
    pred$mean <- as.numeric(pred$mean) + forecast_adjustment
    if (!is.null(pred$lower)) pred$lower <- sweep(as.matrix(pred$lower), 1L, forecast_adjustment, "+")
    if (!is.null(pred$upper)) pred$upper <- sweep(as.matrix(pred$upper), 1L, forecast_adjustment, "+")
  }
  out <- bind_rows(
    as_tibble(pred$lower) %>%
      mutate(horizon = 0:(horizons_to_forecast - 1), lower = TRUE) %>%
      pivot_longer(cols = !c(horizon, lower)),
    as_tibble(pred$upper) %>%
      mutate(horizon = 0:(horizons_to_forecast - 1), lower = FALSE) %>%
      pivot_longer(cols = !c(horizon, lower))
  ) %>%
    rename(quantile = name) %>%
    filter(!(quantile == "0%" & lower)) %>%
    mutate(
      location = meta$location,
      location_name = meta$location_name,
      forecast_date = meta$forecast_date,
      reference_date = meta$reference_date,
      target = paste0(horizon + 1, " wk ahead inc flu hosp"),
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

# MMN fit on z = log(y+1)+1; pass PI on z through to pred_to_quant by shifting to log(y+1).
pred_to_quant_mmn_log1p_shift <- function(pred, meta, forecast_adjustment = NULL) {
  if (!is.null(forecast_adjustment)) {
    if (length(forecast_adjustment) != horizons_to_forecast) {
      stop("forecast_adjustment must have length horizons_to_forecast", call. = FALSE)
    }
    pred$mean <- as.numeric(pred$mean) + forecast_adjustment
    if (!is.null(pred$lower)) pred$lower <- sweep(as.matrix(pred$lower), 1L, forecast_adjustment, "+")
    if (!is.null(pred$upper)) pred$upper <- sweep(as.matrix(pred$upper), 1L, forecast_adjustment, "+")
  }
  pred$mean <- pred$mean - 1
  if (!is.null(pred$lower)) pred$lower <- pred$lower - 1
  if (!is.null(pred$upper)) pred$upper <- pred$upper - 1
  pred_to_quant(pred, meta, log_scale = TRUE, forecast_adjustment = NULL)
}

# MMdN on u = exp(deseasoned_log); recompose to log(y+1) before pred_to_quant(..., log_scale=TRUE).
recompose_mmn_on_exp_deseason <- function(pred, adj_log) {
  pred$mean <- log(pmax(as.numeric(pred$mean), 1e-300)) + adj_log
  if (!is.null(pred$lower)) {
    pred$lower <- log(pmax(as.matrix(pred$lower), 1e-300)) + adj_log
  }
  if (!is.null(pred$upper)) {
    pred$upper <- log(pmax(as.matrix(pred$upper), 1e-300)) + adj_log
  }
  pred
}
