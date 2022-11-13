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

formals(easy_edm)$verbosity <- 0
formals(easy_edm)$showProgressBar <- FALSE


test_that("Simple manifolds", {
  start.time <- Sys.time()
  
  obs <- 500
  map <- logistic_map(obs)

  x <- map$x
  y <- map$y

  df <- data.frame(list(x = x, y = y))

  print("Starting test 1")
  # Check that passing the data via a dataframe works.
  xCCMCausesY <- easy_edm("x", "y", data = df)
  testthat::expect_true(xCCMCausesY)
  print(Sys.time()- start.time)
  
  # yCCMCausesX <- easy_edm("y", "x", data = df)
  # testthat::expect_true(yCCMCausesX)

  print("Starting test 2")
  # Check that passing the raw data is also fine.
  xCCMCausesY <- easy_edm(x, y)
  testthat::expect_true(xCCMCausesY)
  print(Sys.time()- start.time)
})


test_that("Chicago dataset", {
  start.time <- Sys.time()
  
  chicagoURL <- url("https://github.com/EDM-Developers/fastEDM/raw/master/vignettes/chicago.csv")
  chicago <- read.csv(chicagoURL)

  print("Starting test 1")
  crimeCCMCausesTemp <- easy_edm("Crime", "Temperature", data = chicago)
  testthat::expect_false(crimeCCMCausesTemp)
  print(Sys.time()- start.time)
  
  print("Starting test 2")
  tempCCMCausesCrime <- easy_edm("Temperature", "Crime", data = chicago)
  testthat::expect_true(tempCCMCausesCrime)
  print(Sys.time()- start.time)
  
  print("Starting test 3")
  # Check that the results still hold up if we don't normalize the inputs
  crimeCCMCausesTemp <- easy_edm("Crime", "Temperature", data = chicago, normalize = FALSE)
  testthat::expect_false(crimeCCMCausesTemp)
  print(Sys.time()- start.time)
  
  print("Starting test 4")
  tempCCMCausesCrime <- easy_edm("Temperature", "Crime", data = chicago, normalize = FALSE)
  testthat::expect_true(tempCCMCausesCrime)
  print(Sys.time()- start.time)
})
