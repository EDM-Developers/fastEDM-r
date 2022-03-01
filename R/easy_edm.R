
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
easy_edm <- function(cause, effect, time=NULL, data=NULL,
                     direction="oneway", verbosity=1, showProgressBar=NULL,
                     normalize=TRUE) {
  
  if (is.null(showProgressBar)) {
    showProgressBar = verbosity > 0
  }
  
  # If doing a rigorous check, begin by seeing if the cause & effect
  # variable appear to be non-linear dynamical system outputs, or just
  # random noise.
  #TODO
  
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
  
  res <- edm(t, y, E=seq(3, 10), verbosity=0, showProgressBar=showProgressBar)
  
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
  
  # Find the maximum library size using S-map and this E selection
  res <- edm(t, y, E=E_best, algorithm="smap", full=TRUE, saveManifolds=TRUE,
             verbosity=0, showProgressBar=showProgressBar)
  libraryMax <- nrow(res$Ms[[1]])
  
  if (verbosity > 0) {
    cli::cli_alert_info("The maximum library size we can use is {libraryMax}.")
  }
  
  # Set up a grid of library sizes to run the cross-mapping over.
  if (libraryMax >= 500) {
    libraryStart = 100
  } else {
    libraryStart = 10
  }
  libraries <- ceiling(seq(libraryStart, libraryMax, length.out=25))
  
  # Next run the convergent cross-mapping (CCM), using the effect to predict the cause.
  res <- edm(t, y, x, E=E_best, library=libraries, algorithm="smap", k=Inf, shuffle=TRUE,
             verbosity=0, showProgressBar=showProgressBar)
      
  # Make some rough guesses for the Monster exponential fit coefficients
  ccmRes <- res$summary
  firstLibrary <- utils::head(ccmRes$library, 1)
  firstRho <- utils::head(ccmRes$rho, 1)
  finalRho <- utils::tail(ccmRes$rho, 1)

  gammaGuess <- 0.001
  rhoInfinityGuess <- finalRho
  alphaGuess <- (firstRho - rhoInfinityGuess) / exp(-gammaGuess*firstLibrary)
  
  monsterFitStart <- list(
    alpha=alphaGuess, gamma=gammaGuess, rhoInfinity=rhoInfinityGuess
  )
  
  monsterFit <- tryCatch(
    expr = {
      monsterFit <- stats::nls(rho ~ alpha*exp(-gamma*library) + rhoInfinity,
                               data=ccmRes, start=monsterFitStart)
      list(
        alpha=stats::coef(monsterFit)[[1]],
        gamma=stats::coef(monsterFit)[[2]],
        rhoInfinity=stats::coef(monsterFit)[[3]]
      )
    },
    error=function(cond) {
      if (verbosity > 0) {
        cli::cli_alert_danger("Couldn't fit an exponential curve to the rho-L values.")
        cli::cli_alert_danger("This may be a sign that these EDM results are not very reliable.")
      }
      return(list(alpha=NA, gamma=NA, rhoInfinity=finalRho))
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
  
  
  if (givenTimeSeriesNames) {
    alert("{causalSummary} of CCM causation from {cause} to {effect} found.")
  } else {
    alert("{causalSummary} of CCM causation found.")
  }
 
  return(monsterFit$rhoInfinity > 0.5) 
}