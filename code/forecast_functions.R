require(tidyverse)

# fit_nbinom <- function(x, if_all_zeros=NA) {
#   if (all(x == 0)) {
#     if (is.na(if_all_zeros)) return(NA)
#     else if (if_all_zeros == 'add1') x <- c(x, 1)
#   }
#   start_size <- min(abs(mean(x)^2 / (var(x) - mean(x))), 1e4) # Poisson-like if negative or high number
#   this_fit <- fitdist(x, 'nbinom', method='mle', start = list(mu=median(x), size = start_size))$estimate
#   tibble(mu=this_fit[['mu']], size=this_fit[['size']])
# }

outlier_check <- function(mu, size, x, weights, if_all_zeros, start=T, threshold = 0.01) {
  outlier_detected <- F
  if (start) {
    x <- x[1:ceiling(length(x)/2)]
    weights <- weights[1:ceiling(length(weights)/2)]
  } else {
    x <- x[floor(length(x)/2):length(x)]
    weights <- weights[floor(length(weights)/2):length(weights)]
  }
  if (all(x == 0) & if_all_zeros == 'add1') {
    x <- c(x, 1)
    weights <- c(weights, 1)
  }
  mu_check <- sum(x * weights) / sum(weights)
  varx <- max(mu_check, var(x))
  size_check <- mu_check^2 / (varx - mu_check)
  p_check <- pnbinom(mu, mu=mu_check, size=size_check)
  if (p_check < threshold | (1-threshold) < p_check) {
    return(tibble(mu = mu_check, size = size_check, outlier_detected = T))
  } else {
    return(tibble(mu = mu, size = size, outlier_detected = F))
  }
}

fit_nbinom <- function(x, if_all_zeros=NA, weights=NA, outlier_check=F) {
  if (all(is.na(weights))) weights <- rep(1, length(x))
  if (all(x == 0)) {
    if (is.na(if_all_zeros)) return(NA)
    else if (if_all_zeros == 'add1') {
      x <- c(x, 1)
      weights <- c(weights, 1)
    }
  }
  #print(paste(length(x), length(weights)))
  if (length(x) != length(weights)) stop('x and weights have different lengths')
  mu <- sum(x * weights) / sum(weights)
  varx <- max(mu, var(x))
  size <- mu^2 / (varx - mu)
  if (outlier_check) { # check for divergence in the later time points (e.g., week 3 if weekly)
    check_start <- outlier_check(mu=mu, size=size, x=x, weights=weights, if_all_zeros=if_all_zeros)
    check_end <- outlier_check(mu=mu, size=size, x=x, weights=weights, if_all_zeros=if_all_zeros, start=F)
    if (check_start$outlier_detected & check_end$outlier_detected) {
      stop('check_start & check_end')
    } else if (check_start$outlier_detected) {
      return(check_start)
    } else return(check_end)
  }
  else return(tibble(mu=mu, size=size))
}

nbinom_to_quant <- function(mu, size, quantiles) {
  tibble(
    quantile = quantiles, 
    value = qnbinom(p=quantiles, mu=mu, size=size))
}

# normal on x or log(x + 1)
norm_to_quant <- function(mean, sd, quantiles, log_transformed=T) {
  tibble(
    quantile = quantiles, 
     value = qnorm(p=quantiles, mean=mean, sd=sd)) %>%
     mutate(
       value = if (log_transformed) exp(value) - 1 else value,
       value = ifelse(value > mean, ceiling(value), floor(value)),
       value = pmax(0, value))
}

# samples from differences
dsample_to_quant <- function(mean, diffs, quantiles, log_transformed=F) {
  diffs <- c(diffs, -diffs) # symmetrize
  diff_quants <- quantile(diffs, probs = quantiles)
  values <- mean + diff_quants
  if (log_transformed) {
    values <- exp(values) - 1
  }
  tibble(
      quantile = quantiles, 
      value = values) %>%
    mutate(
      value = pmax(value, 0), # force non-negative incidence
      value = ifelse(value > mean, ceiling(value), floor(value)))
}

dsample_to_bin <- function(mean, diffs, bins, n_sample=2*length(diffs)) {
  if (min(bins$bin_start_incl) != 0 | max(bins$bin_end_notincl) != 1) {
    stop('min(bin_start_incl) must be zero and max(bin_end_notincl) must be one')
  }
  diffs <- c(diffs, -diffs) # symmetrize
  sampled_values <- pmin(pmax(mean + diffs, 0), 1) # force to 0-1
  sample_ecdf <- ecdf(sampled_values)
  tibble(
    bin_start_incl = bins$bin_start_incl,
    bin_end_notincl = bins$bin_end_notincl,
    value = sample_ecdf(bins$bin_end_notincl) - sample_ecdf(bins$bin_start_incl)
  ) %>%
  mutate(value = ifelse(bins$bin_start_incl == 0, value + sample_ecdf(0), value))
}

# logit
fit_logitnorm <- function(x, zero_val=NA, weights=rep(1, length(x))) {
  if (all(x == 0)) stop('All values are zero')
  if (any(x == 0) & is.na(zero_val)) { 
    stop('At least one 0 detect (logit(0) = -Inf), set zero_val')
  } else {
    x[x == 0] <- zero_val
  }
  if (length(x) != length(weights)) stop('x and weights have different lengths')
  #if (any(x == 0)) x[x == 0] <- min(x[x != 0])
  x_logit <- boot::logit(x)
  mu <- sum(x_logit * weights) / sum(weights)
  tibble(logit_mean=mu, logit_sd=sd(x_logit))
}

logitnorm_to_bin <- function(logit_mean, logit_sd, bins) {
  if (min(bins$bin_start_incl) != 0 | max(bins$bin_end_notincl) != 1) {
    stop('min(bin_start_incl) must be zero and max(bin_end_notincl) must be one')
  }
  tibble(
    bin_start_incl = bins$bin_start_incl,
    bin_end_notincl = bins$bin_end_notincl,
    value = pnorm(boot::logit(bins$bin_end_notincl), mean=logit_mean, sd=logit_sd) - 
      pnorm(boot::logit(bins$bin_start_incl), mean=logit_mean, sd=logit_sd)
  )
}

logitnorm_to_quant <- function(logit_mean, logit_sd, quantiles) {
  tibble(
    quantile = quantiles, 
    value = qnorm(p=quantiles, mean=logit_mean, sd=logit_sd)) %>%
    mutate(value = boot::inv.logit(value))
}
