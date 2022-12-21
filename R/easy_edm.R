
#' easy_edm
#' 
#' This is an automated workflow for performing causal analysis on the
#' supplied time series using empirical dynamical modelling (EDM) techniques.
#' It is intended to hide all the common steps of an EDM analysis, and should
#' work on most datasets.
#' 
#' It may be the case that your data requires a custom analysis, so this
#' function can be used as a helpful starting point from which to create a
#' specialised analysis using the `edm` function directly.
#' 
#' Warning: While the `edm` functionality is well-tested and ready for use,
#' this `easy_edm` automated analysis is still a work-in-progress.
#'
#' @param cause The causal time series (as a string or a vector).
#'
#' @param effect The effect time series (as a string or a vector).
#'
#' @param time For non-regularly sampled time series, the sampling times
#' must be supplied here (as a string or a vector).
#'
#' @param data If a dataframe is supplied here, then the cause, effect & time
#' arguments must be the column names of the relevant time series as strings.
#'
#' @param direction A string specifying whether we are checking a one
#' directional causal effect or whether to test the reverse direction at the
#' same time (work in progress!).
#'
#' @param verbosity The level of detail in the output.
#'
#' @param showProgressBar Whether or not to print out a progress bar during the computations.
#'
#' @param normalize Whether to normalize the inputs before starting EDM.
#'
#' @returns A Boolean indicating that evidence of causation was found.
#' @export
#' @example man/examples/logistic-map-example.R
#
easy_edm <- function(cause, effect, time = NULL, data = NULL,
                     direction = "oneway", verbosity = 1, showProgressBar = NULL,
                     normalize = TRUE) {

  # !! Parameterise these values later
  convergence_method <- "quantile"
  max_theta <- 5
  num_thetas <- 100
  theta_reps <- 20
  max_lag <- 5

  if (is.null(showProgressBar)) {
    showProgressBar <- verbosity > 0
  }

  # Convert time series to vectors (they can be supplied as columns of a dataframe).
  inputs <- preprocess_inputs(data, cause, effect, time, verbosity, normalize)
  t <- inputs$t
  x <- inputs$x
  y <- inputs$y

  # Find optimal E (embedding dimension) of the causal variable using simplex projection
  E_best <- find_embedding_dimension(t, x, verbosity, showProgressBar)

  # Test for non-linearity using S-Map
  test_output <- test_nonlinearity(t, x, E_best, max_theta, num_thetas, theta_reps,
                                  verbosity, showProgressBar)
  optTheta    <- test_output$theta
  isNonLinear <- test_output$outcome

  # Lags the y (effect) time series by the optimal value or differences the series if it was linear
  yOpt <- get_optimal_effect(t, x, y, E_best, verbosity, showProgressBar, isNonLinear, optTheta, max_lag)

  # Get max library size
  libraryMax <- get_max_library(t, x, yOpt, E_best, verbosity, showProgressBar)

  # Test for causality using CCM
  if (convergence_method == "parametric") {
      result <- test_convergence_monster(t, x, yOpt, E_best, libraryMax, optTheta, verbosity, showProgressBar)
  }
  else if (convergence_method == "quantile") {
      result <- test_convergence_dist(t, x, yOpt, E_best, libraryMax, optTheta, verbosity, showProgressBar)
  }
  else {
      result <- test_convergence_dist(t, x, yOpt, E_best, libraryMax, optTheta, verbosity, showProgressBar)
  }

  foundEvidence <- result != "No evidence"
  if (foundEvidence) {
    alert <- cli::cli_alert_success
  } else {
    alert <- cli::cli_alert_danger
  }

  givenTimeSeriesNames <- !is.null(data)
  if (givenTimeSeriesNames) {
    alert("{result} of CCM causation from {cause} to {effect} found.")
  } else {
    alert("{result} of CCM causation found.")
  }

  return(foundEvidence)
}

# ---------------------------------------------------------------------------------------
preprocess_inputs <- function(data, cause, effect, time, verbosity, normalize) {
  # First find out the embedding dimension of the causal variable
  givenTimeSeriesNames <- !is.null(data)
  if (givenTimeSeriesNames) {
    if (verbosity > 0) {
      cli::cli_alert_info("Pulling the time series from the supplied dataframe.")
    }

    if (!(cause %in% colnames(data))) {
      cli::cli_alert_danger("{cause} is not a column in the supplied dataframe.")
      return(1)
    }
    if (!(effect %in% colnames(data))) {
      cli::cli_alert_danger("{effect} is not a column in the supplied dataframe.")
      return(1)
    }
    if (!is.null(time) && !(time %in% colnames(data))) {
      cli::cli_alert_danger("{time} is not a column in the supplied dataframe.")
      return(1)
    }

    x <- data[[cause]]
    y <- data[[effect]]
    t <- if (is.null(time)) seq(length(x)) else data[[time]]
  } else {
    if (verbosity > 0) {
      cli::cli_alert_info("Using supplied time series vectors.")
    }
    x <- cause
    y <- effect
    t <- if (is.null(time)) seq(length(x)) else time
  }

  if (length(t) != length(x) || length(t) != length(y)) {
    cli::cli_alert_danger("Time series are not the same length.")
    return(1)
  }

  if (verbosity > 0) {
    cli::cli_alert_info("Number of observations is {length(t)}")
  }

  if (normalize) {
    if (verbosity > 0) {
      cli::cli_alert_info("Normalizing the supplied time series")
    }
    x <- scale(x)
    y <- scale(y)
  }

  return(data.frame(t=t, x=x, y=y))
}

# ---------------------------------------------------------------------------------------
# Find optimal E (embedding dimension) of the causal variable using simplex projection
find_embedding_dimension <- function(t, x, verbosity, showProgressBar) {

  res <- edm(t, x, E = seq(3, 10), verbosity = 0, showProgressBar = showProgressBar)

  if (res$rc > 0) {
    cli::cli_alert_danger("Search for optimal embedding dimension failed.")
    return(2)
  }

  if (is.null(res$summary$rho) | length(res$summary$rho) == 0) {
    cli::cli_alert_danger("Search for optimal embedding dimension failed (2).")
    return(3)
  }

  E_best <- res$summary$E[which.max(res$summary$rho)]

  if (verbosity > 0) {
    cli::cli_alert_success("Found optimal embedding dimension E to be {E_best}.")
  }

  return(E_best)
}

# ---------------------------------------------------------------------------------------
# Test for non-linearity using S-Map
test_nonlinearity <- function(t, x, E_best, max_theta, num_thetas, theta_reps, verbosity, showProgressBar) {

  theta_values <- seq(0, max_theta, length.out=1+num_thetas)

  res <- edm(t, x, E = E_best, theta = theta_values, algorithm="smap", k=Inf,
             verbosity = 0, showProgressBar = showProgressBar)

  optIndex <- which(res$summary$rho==max(res$summary$rho))
  optRho   <- res$summary$rho[optIndex]
  optTheta <- res$summary$theta[optIndex]

  if (verbosity > 0) {
    cli::cli_alert_success("Found optimal theta to be {optTheta}, with rho = {optRho}.")
  }

  resBase <- edm(t, x, E = E_best, theta = 0, verbosity = 0, numReps = theta_reps,
                 k=20, algorithm = "smap", showProgressBar = showProgressBar)
  resOpt  <- edm(t, x, E = E_best, theta = optTheta, verbosity = 0, numReps = theta_reps,
                 k=20, algorithm = "smap", showProgressBar = showProgressBar)

  sampleBase <- resBase$stats$rho
  sampleOpt  <- resOpt$stats$rho

  ksOut  <- stats::ks.test(sampleOpt, sampleBase, alternative="less")
  ksStat <- ksOut$statistic
  ksPVal <- ksOut$p.value

  if (verbosity > 0) {
    cli::cli_alert_success("Found Kolmogorov-Smirnov test statistic to be {ksStat} with p-value={ksPVal}.")
  }

  isNonLinear <- ksPVal < 0.05

  return(data.frame(theta=optTheta, outcome=isNonLinear))
}

# ---------------------------------------------------------------------------------------
tslag <- function(t, x, lag = 1, dt = 1) {
  l.x <- rep(NA, length(t))
  for (i in seq_along(t)) {
    lagged_t <- t[i] - lag * dt
    if (!is.na(lagged_t) && lagged_t %in% t) {
      l.x[i] <- x[which(t == lagged_t)]
    }
  }
  return(l.x)
}

# ---------------------------------------------------------------------------------------
# Find optimal lag for the y (effect) time series
get_optimal_effect <- function(t, x, y, E_best, verbosity, showProgressBar, isNonLinear, theta, maxLag) {
    lagRhos = data.frame(lag=double(), rho=double())
    for (i in seq(-maxLag, maxLag + 1, 1)) {
      res = edm(t, x, tslag(t, y, i), E = E_best, theta = theta,
                algorithm="simplex", k=Inf, verbosity = 0, showProgressBar = showProgressBar)
      lagRhos[nrow(lagRhos) + 1,] = c(i, res$summary$rho)
    }

    # Sort by Rho values
    lagRhos <- lagRhos[order(-lagRhos$rho),]

    optLag <- lagRhos$lag[1]
    optRho <- round(lagRhos$rho[1], 5)

    if (verbosity > 0)
      cli::cli_alert_info('Found optimal lag to be {optLag} with rho={optRho}')

    # If retro-causality is spotted,
    # default to best positive lag and print warning
    invalidLag = optLag < 0
    if (invalidLag) {
      validRhos <- lagRhos[lagRhos >= 0,]
      optLag <- validRhos$lag[1]
      optRho <- round(validRhos$rho[1], 5)

      if (verbosity > 0)
        cli::cli_alert_info(
          'Potential retrocausality, using lag of {optLag} with rho={optRho}'
        )
    }

    yLag <- tslag(t, y, optLag)
    if (isNonLinear) {
      # Transform y and data to match optimal lagged series
      yOpt <- yLag

      if (verbosity > 0) {
        cli::cli_alert_info("Lagging time series using optimal lag of {optLag}")
      }
    }
    else {
      # Difference y and data by optimal lagged series
      yOpt <- y - yLag

      if (verbosity > 0)
        cli::cli_alert_info('Differencing time series due to failed nonlinearity test (lag={optLag})')
    }

    return(yOpt)
}

# ---------------------------------------------------------------------------------------
# Find maximum library size
get_max_library <- function(t, x, y, E_best, verbosity, showProgressBar) {
  res <- edm(t, y,
             E = E_best, algorithm = "simplex", full = TRUE, saveManifolds = TRUE,
             verbosity = 0, showProgressBar = showProgressBar
  )
  libraryMax <- nrow(res$Ms[[1]])
  return(libraryMax)
}

# ---------------------------------------------------------------------------------------
# Test for causality using CCM
cross_mapping <- function(t, x, y, E_best, verbosity, showProgressBar) {
  # Find the maximum library size using S-map and this E selection
  res <- edm(t, y,
    E = E_best, algorithm = "simplex", full = TRUE, saveManifolds = TRUE,
    verbosity = 0, showProgressBar = showProgressBar
  )
  libraryMax <- nrow(res$Ms[[1]])

  if (verbosity > 0) {
    cli::cli_alert_info("The maximum library size we can use is {libraryMax}.")
  }

  # Set up a grid of library sizes to run the cross-mapping over.
  if (libraryMax >= 500) {
    libraryStart <- 100
  } else {
    libraryStart <- 10
  }
  libraries <- ceiling(seq(libraryStart, libraryMax, length.out = 25))

  # Next run the convergent cross-mapping (CCM), using the effect to predict the cause.
  res <- edm(t, y, x,
    E = E_best, library = libraries, algorithm = "simplex", k = Inf, shuffle = TRUE,
    verbosity = 0, showProgressBar = showProgressBar
  )

  return(res)
}

# ---------------------------------------------------------------------------------------
# Test for convergence using parametric test (Monster)
test_convergence_monster <- function(t, x, y, E_best, libraryMax, theta, verbosity, showProgressBar) {
  # Perform cross-mapping (CCM)
  ccm <- cross_mapping(t, x, y, E_best, verbosity, showProgressBar)

  # Make some rough guesses for the Monster exponential fit coefficients
  ccmRes <- ccm$summary
  firstLibrary <- utils::head(ccmRes$library, 1)
  firstRho <- utils::head(ccmRes$rho, 1)
  finalRho <- utils::tail(ccmRes$rho, 1)

  gammaGuess <- 0.001
  rhoInfinityGuess <- finalRho
  alphaGuess <- (firstRho - rhoInfinityGuess) / exp(-gammaGuess * firstLibrary)

  monsterFitStart <- list(
    alpha = alphaGuess, gamma = gammaGuess, rhoInfinity = rhoInfinityGuess
  )

  monsterFit <- tryCatch(
    expr = {
      monsterFit <- stats::nls(rho ~ alpha * exp(-gamma * library) + rhoInfinity,
        data = ccmRes, start = monsterFitStart
      )
      list(
        alpha = stats::coef(monsterFit)[[1]],
        gamma = stats::coef(monsterFit)[[2]],
        rhoInfinity = stats::coef(monsterFit)[[3]]
      )
    },
    error = function(cond) {
      if (verbosity > 0) {
        cli::cli_alert_danger("Couldn't fit an exponential curve to the rho-L values.")
        cli::cli_alert_danger("This may be a sign that these EDM results are not very reliable.")
      }
      return(list(alpha = NA, gamma = NA, rhoInfinity = finalRho))
    }
  )

  if (verbosity > 1) {
    cli::cli_alert_info("The CCM fit is (alpha, gamma, rhoInfinity) = ({signif(monsterFit$alpha, 2)}, {signif(monsterFit$gamma, 2)}, {signif(monsterFit$rhoInfinity, 2)}).")
  } else if (verbosity == 1) {
    cli::cli_alert_info("The CCM final rho was {signif(monsterFit$rhoInfinity, 2)}")
  }

  if (monsterFit$rhoInfinity > 0.7) {
    causalSummary <- "Strong evidence"
  } else if (monsterFit$rhoInfinity > 0.5) {
    causalSummary <- "Some evidence"
  } else {
    causalSummary <- "No evidence"
  }

  return(causalSummary)
}

# ---------------------------------------------------------------------------------------
# Test for convergence by comparing the distribution of rho at a small library size and
# a sampled rho at the maximum library size.
test_convergence_dist <- function(t, x, y, E_best, libraryMax, theta, verbosity, showProgressBar, numReps = 1000) {

  librarySmall <- max(E_best + 2, floor(libraryMax/10))

  distRes <- edm(t, y, x, E = E_best, library = librarySmall, numReps = numReps,
                 theta = theta, algorithm = "simplex", k = Inf, shuffle = TRUE,
                 verbosity = 0, showProgressBar = showProgressBar)
  dist <- distRes$stats$rho

  finalRes <- edm(t, y, x, E = E_best, library = libraryMax, theta = theta,
                  algorithm = "simplex", k = Inf, shuffle = TRUE,
                  verbosity = 0, showProgressBar = showProgressBar)
  finalRho <- finalRes$summary$rho

  quantiles <- stats::quantile(dist, c(0.95, 0.975))

  q975 <- quantiles[[1]]
  q95  <- quantiles[[2]]

  rhoQuantile = length(dist[dist < finalRho]) / length(dist)

  if (verbosity >= 1 ){
    cli::cli_alert_info("At library=E+2, found rho quantiles of {round(q975, 5)} (0.975) and {round(q95, 5)} (0.95)")
    cli::cli_alert_info("At library=max, found final rho was {round(finalRho, 5)}, i.e. quantile={rhoQuantile}")
  }

  if (finalRho > q975) {
    causalSummary = "Strong evidence"
  }
  else if (finalRho > q95) {
    causalSummary = "Some evidence"
  }
  else {
    causalSummary = "No evidence"
  }

  return(causalSummary)
}
