
#' easy_edm
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
#' @example man/chicago-easy-edm-example.R
#
easy_edm <- function(cause, effect, time = NULL, data = NULL,
                     direction = "oneway", verbosity = 1, showProgressBar = NULL,
                     normalize = TRUE) {

  # !! Parameterise these values later
  max_theta <- 5
  num_thetas <- 100
  theta_reps <- 20
  convergence_method <- "parametric"

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
  optTheta <- test_nonlinearity(t, x, E_best, max_theta, num_thetas, theta_reps, verbosity, showProgressBar)

  # Perform cross-mapping (CCM)
  res <- cross_mapping(t, x, y, E_best, verbosity, showProgressBar)

  # Test for causality using CCM
  if (convergence_method == "parametric") {
    conv_test <- test_convergence_monster
  }
  else if (convergence_method == "hypothesis") {
      conv_test <- test_convergence_monster # Replace this later
  }
  else {
      conv_test <- test_convergence_monster # Replace this later
  }

  outcome <- conv_test(res, data, cause, effect, verbosity)

  return(outcome)
}

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
  debug = FALSE
  
  max_theta <- 5; theta_step <- 0.05; theta_reps <- 20;

  theta_values <- seq(0, max_theta, theta_step)

  res <- edm(t, x, E = E_best, theta = theta_values, algorithm="smap", k=Inf,
             verbosity = 0, showProgressBar = showProgressBar)

  optIndex <- which(res$summary$rho==max(res$summary$rho))
  optRho   <- res$summary$rho[optIndex]
  optTheta <- res$summary$theta[optIndex]

  if (verbosity > 0 || debug) {
    cli::cli_alert_success("Found optimal theta to be {optTheta}, with rho = {optRho}.")
  }

  resBase <- edm(t, x, E = E_best, theta = 0, verbosity = 0, numReps = theta_reps,
                 k=20, algorithm = "smap", showProgressBar = showProgressBar)
  resOpt  <- edm(t, x, E = E_best, theta = optTheta, verbosity = 0, numReps = theta_reps,
                 k=20, algorithm = "smap", showProgressBar = showProgressBar)

  sampleBase <- resBase$stats$rho
  sampleOpt  <- resOpt$stats$rho

  ksOut  <- ks.test(sampleOpt, sampleBase, alternative="less")
  ksStat <- ksOut$statistic
  ksPVal <- ksOut$p.value

  if (verbosity > 0 || debug) {
    cli::cli_alert_success("Found Kolmogorov-Smirnov test statistic to be {ksStat} with p-value={ksPVal}.")
  }

  return(optTheta)
}

# ---------------------------------------------------------------------------------------
# Test for causality using CCM
cross_mapping <- function(t, x, y, E_best, verbosity, showProgressBar) {
  # Find the maximum library size using S-map and this E selection
  res <- edm(t, y,
    E = E_best, algorithm = "smap", full = TRUE, saveManifolds = TRUE,
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
    E = E_best, library = libraries, algorithm = "smap", k = Inf, shuffle = TRUE,
    verbosity = 0, showProgressBar = showProgressBar
  )

  return(res)
}

# ---------------------------------------------------------------------------------------
# Test for convergence using parametric test (Monster)
test_convergence_monster <- function(res, data, cause, effect, verbosity) {
  # Make some rough guesses for the Monster exponential fit coefficients
  ccmRes <- res$summary
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
    alert <- cli::cli_alert_success
  } else if (monsterFit$rhoInfinity > 0.5) {
    causalSummary <- "Some evidence"
    alert <- cli::cli_alert_success
  } else {
    causalSummary <- "No evidence"
    alert <- cli::cli_alert_danger
  }

  givenTimeSeriesNames <- !is.null(data)
  if (givenTimeSeriesNames) {
    alert("{causalSummary} of CCM causation from {cause} to {effect} found.")
  } else {
    alert("{causalSummary} of CCM causation found.")
  }

  return(monsterFit$rhoInfinity > 0.5)
}
