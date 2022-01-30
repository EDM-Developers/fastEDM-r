
#' easy_edm
#'
#' @param cause The first time series in the causal analysis
#' 
#' @param effect The second time series in the causal analysis
#' 
#' @param time If time is not uniformly sampled, then it must be supplied here.
#' 
#' @param data If a dataframe is supplied here, then cause, effect & time must
#' be strings containing the column names of the relevant time series.
#' 
#' @param direction A string specifying whether we are checking a one
#' directional causal effect or whether to test the reverse direction at the
#' same time.
#' 
#' @param verbosity The level of detail in the output.
#' 
#' @param normalize Whether to normalize the inputs before starting EDM.
#' 
#' @returns An integer error/return code (success is 0)
#' @export
#' @examples
#'  library(fastEDM)
#'  library(readr)
#'  data <- url("https://raw.githubusercontent.com/EDM-Developers/EDM/master/test/chicago.csv")
#'  
#'  chicago <- read_csv(data, col_types = cols(crime = col_double()))
#'  chicago <- head(chicago, 500) # Just to speed up the example
#'  easy_edm("crime", "temp", data=chicago)
#
easy_edm <- function(cause, effect, time=NULL, data=NULL,
                     direction="oneway", verbosity=1, normalize=TRUE) {
  
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
    x <- scale(x)
    y <- scale(y)
  }
  
  res <- edm(t, y, E=seq(3, 10), verbosity=verbosity)
  
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
  
  # Find the maximum library size using this E selection
  res <- edm(t, y, E=E_best, full=TRUE, saveManifolds=TRUE, verbosity=verbosity)
  libraryMax <- length(res$Ms[[1]])
  
  if (verbosity > 0) {
    cli::cli_alert_info("The maximum library size we can use is {libraryMax}.")
  }
  
  # Next do causal cross-mapping (CCM) from the cause to the effect
  libraries <- ceiling(seq(10, libraryMax, length.out=25))
  
  res <- edm(t, y, x, E=E_best, library=libraries, verbosity=verbosity)
      
  # Make some rough guesses for the Monster exponential fit coefficients
  ccmRes <- res$summary
  firstLibrary <- head(ccmRes$library, 1)
  firstRho <- head(ccmRes$rho, 1)
  finalRho <- tail(ccmRes$rho, 1)

  gammaGuess <- 0.01
  rhoInfinityGuess <- finalRho
  alphaGuess <- (firstRho - rhoInfinityGuess) / exp(-gammaGuess*firstLibrary)
  
  monsterFit <- nls(rho ~ alpha*exp(-gamma*library) + rhoInfinity, data=ccmRes,
                    start = list(alpha=alphaGuess, gamma=gammaGuess, rhoInfinity=rhoInfinityGuess))

  alphaFitted <- coef(monsterFit)[[1]]
  gammaFitted <- coef(monsterFit)[[2]]
  rhoInfinityFitted <- coef(monsterFit)[[3]]
  
  if (verbosity > 0) {
    cli::cli_alert_info("The CCM fit is (alpha, gamma, rhoInfinity) = ({signif(alphaFitted, 2)}, {signif(gammaFitted, 2)}, {signif(rhoInfinityFitted, 2)}).")
  }

  
  if (rhoInfinityFitted > 0.7) {
    if (givenTimeSeriesNames) {
      cli::cli_alert_success("Strong evidence that {cause} causes {effect}.")
    } else {
      cli::cli_alert_success("Strong evidence of causal link.")
    }
  } else {
    if (givenTimeSeriesNames) {
      cli::cli_alert_danger("No causal link from {cause} to {effect} found.")
    } else {
      cli::cli_alert_danger("No causal link found.")
    }
  }
  
  return(rhoInfinityFitted > 0.7) 
}