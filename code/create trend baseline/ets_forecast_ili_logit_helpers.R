#!/usr/bin/env Rscript
# ETS helpers for ILI using logit(weighted_ili / 100) transform.

source("code/ili_forecast_common.R")

horizons_to_forecast <- 12L
quantiles_ili <- c(
  0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
  0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99
)
forecast_levels <- 100 * c(1 - 2 * quantiles_ili[1:12])
STL_PERIOD <- 52L
LOGIT_EPS <- 1e-6

ili_to_prop <- function(y_pct) {
  pmin(pmax(y_pct / 100, LOGIT_EPS), 1 - LOGIT_EPS)
}

stl_future_seasonal <- function(seas_vec, period, h_max) {
  last_cycle <- tail(seas_vec, period)
  idx <- ((seq_len(h_max) - 1L) %% period) + 1L
  as.numeric(last_cycle[idx])
}

stl_prefilter_logit <- function(y_pct, period, h_max) {
  n <- length(y_pct)
  if (n < 2L * period) return(NULL)
  x <- boot::logit(ili_to_prop(y_pct))
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

pred_to_bin_logit <- function(pred, meta, bins, forecast_adjustment = NULL) {
  if (!is.null(forecast_adjustment)) {
    if (length(forecast_adjustment) != horizons_to_forecast)
      stop("forecast_adjustment must have length horizons_to_forecast", call. = FALSE)
    pred$mean  <- as.numeric(pred$mean) + forecast_adjustment
    if (!is.null(pred$upper))
      pred$upper <- sweep(as.matrix(pred$upper), 1L, forecast_adjustment, "+")
  }
  means    <- as.numeric(pred$mean)
  upper_98 <- as.numeric(pred$upper[, "98%"])
  sigmas   <- (upper_98 - means) / qnorm(0.99)

  out <- tibble()
  for (h in seq_len(horizons_to_forecast)) {
    bin_rows <- logitnorm_to_bin(means[h], sigmas[h], bins) %>%
      mutate(
        location        = meta$location,
        forecast_date   = meta$forecast_date,
        target          = paste0(h, " wk ahead"),
        target_end_date = meta$target_end_date_base + (h - 1L) * 7L
      )
    out <- bind_rows(out, bin_rows)
  }
  out
}

pred_to_quant_logit <- function(pred, meta, forecast_adjustment = NULL) {
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
      target = paste0(horizon + 1, " wk ahead"),
      target_end_date = meta$target_end_date_base + horizon * 7,
      quantile = as.numeric(str_remove(quantile, "%")) / 100,
      quantile = ifelse(lower, (1 - quantile) / 2, 1 - (1 - quantile) / 2)
    ) %>%
    mutate(
      value = 100 * boot::inv.logit(value),
      value = pmin(pmax(value, 0), 100)
    )

  out %>%
    mutate(value = round(value, 4)) %>%
    arrange(horizon, quantile) %>%
    dplyr::select(-lower, -horizon)
}
