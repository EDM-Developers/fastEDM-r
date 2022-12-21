
#' logistic_map
#' 
#' This function generates a synthetic time series called the logistic map.
#' It can be used to perform basic tests of EDM causality.
#'
#' @param obs The number of observations to generate.
#'
#' @returns A list containing the synthetic x & y time series.
#' @export
#
logistic_map <- function(obs) {
  r_x <- 3.625
  r_y <- 3.77
  beta_xy <- 0.05
  beta_yx <- 0.4
  tau <- 1
  
  x <- rep(NA, obs)
  y <- rep(NA, obs)
  
  x[1] <- 0.2
  y[1] <- 0.4
  
  for (i in 2:obs) {
    x[i] <- x[i - 1] * (r_x * (1 - x[i - 1]) - beta_xy * y[i - 1])
    y[i] <- y[i - 1] * (r_y * (1 - y[i - 1]) - beta_yx * x[i - tau])
  }
  
  return(list(x = x, y = y))
}
