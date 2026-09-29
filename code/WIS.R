weighted_interval_score <- function(quantile, value, obs_value) {
  if (all(is.na(obs_value))) return(NA)
  actual_value <- obs_value[[1L]]

  value <- value[!is.na(quantile)]
  quantile <- quantile[!is.na(quantile)]

  2 * mean(pmax(
    quantile * (actual_value - value),
    (1 - quantile) * (value - actual_value), 
    na.rm = TRUE))
}

weighted_interval_score_decomp <- function(quantiles, values, obs_value) {
  if (all(is.na(obs_value))) {
    return(data.frame(wis=NA, dispersion=NA, underprediction=NA, overprediction=NA))
  }
  actual_value <- obs_value[1]
  if (any(obs_value != actual_value)) stop('multiple observation values provided')
  if (length(quantiles) != length(values)) stop('numbers of quantiles and values do not match')
  K <- (length(quantiles) - 1)/2
  
  median_value <- values[quantiles == 0.5]
  
  alpha <- 1 - abs(0.5 - quantiles) * 2
  weight = ifelse(quantiles == 0.5, 0.5 * alpha/2, alpha/2)
  half_width <- ifelse(quantiles != 0.5, abs(values - median_value), 0)
  overprediction <- 2/alpha * ifelse(quantiles <= 0.5 & actual_value <= values, values - actual_value, 0)
  underprediction <- 2/alpha * ifelse(quantiles >= 0.5 & actual_value >= values, actual_value - values, 0)

  dispersion <- 1/(K + 1/2) * sum(weight * half_width)
  underprediction <- 1/(K + 1/2) * sum(weight * underprediction)
  overprediction <- 1/(K + 1/2) * sum(weight * overprediction)
  wis <- dispersion + overprediction + underprediction
  
  return(data.frame(wis=wis, dispersion=dispersion, underprediction=underprediction, overprediction=overprediction))
}


