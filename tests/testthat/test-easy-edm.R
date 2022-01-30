
library(readr)

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
    x[i] <- x[i-1] * (r_x * (1 - x[i-1]) - beta_xy * y[i-1])
    y[i] <- y[i-1] * (r_y * (1 - y[i-1]) - beta_yx * x[i-tau])
  }

  return(list(x=x, y=y))
}


formals(edm)$verbosity <- 0

test_that("Simple manifolds", {
  obs <- 500
  map <- logistic_map(obs)
  
  x <- map$x
  y <- map$y

  df <- data.frame(list(x = x, y = y))
  
  # Check that passing the data via a dataframe works.
  xCCMCausesY <- easy_edm("x", "y", data=df)
  testthat::expect_true(xCCMCausesY)
  
  yCCMCausesX <- easy_edm("y", "x", data=df)
  testthat::expect_true(yCCMCausesX)

  # Check that passing the raw data is also fine.
  xCCMCausesY <- easy_edm(x, y)
  testthat::expect_true(xCCMCausesY)
  
})

test_that("Chicago dataset", {
  data <- url("https://raw.githubusercontent.com/EDM-Developers/EDM/master/test/chicago.csv")
  
  chicago <- read_csv(data, col_types = cols(crime = col_double()))
  chicago <- head(chicago, 500) # Just to speed up the example
  
  crimeCCMCausesTemp <- easy_edm("crime", "temp", data=chicago)
  testthat::expect_false(crimeCCMCausesTemp)
  
  tempCCMCausesCrime <- easy_edm("temp", "crime", data=chicago)
  testthat::expect_true(tempCCMCausesCrime)
  
  # Check that the results still hold up if we don't normalize the inputs
  crimeCCMCausesTemp <- easy_edm("crime", "temp", data=chicago, normalize=FALSE)
  testthat::expect_false(crimeCCMCausesTemp)
  
  tempCCMCausesCrime <- easy_edm("temp", "crime", data=chicago, normalize=FALSE)
  testthat::expect_true(tempCCMCausesCrime)
})
