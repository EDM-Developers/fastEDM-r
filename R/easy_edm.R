
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
#' @returns list
#' @export
easy_edm <- function(cause, effect, time=NULL, data=NULL,
                     direction="oneway", verbosity=0) {
  
  # If doing a rigorous check, begin by seeing if the cause & effect
  # variable appear to be non-linear dynamical system outputs, or just
  # random noise.
  #TODO
  
  
  # First find out the embedding dimension of the causal variable
  if (is.null(data)) {
    if (verbosity > 0) {
      cli::cli_alert_info("Using supplied time series vectors.")
    }
    x <- cause
    y <- effect
    t <- if (is.null(time)) seq(length(x)) else time
  } else {
    if (verbosity > 0) {
      cli::cli_alert_info("Pulling the time series from the supplied dataframe.")
    }
    x <- data[[cause]]
    y <- data[[effect]]
    t <- if (is.null(time)) seq(length(x)) else data[[time]]
  }
  
  if (length(t) != length(x) || length(t) != length(y)) {
    cli::cli_alert_danger("Time series are not the same length.")
    return()
  }
  
  if (verbosity > 0) {
    cli::cli_alert_info("Number of observations is {length(t)}")
  }
  
  res <- edm(t, x, E=seq(3, 10), verbosity=verbosity)
  
  if (res$rc > 0) {
    cli::cli_alert_danger("Search for optimal embedding dimension failed.")
    return()
  }
  
  if (is.null(res$summary$rho) | length(res$summary$rho) == 0) {
    cli::cli_alert_danger("Search for optimal embedding dimension failed (2).")
    return()
  }
  
  E_best <- res$summary$E[which.max(res$summary$rho)]
  
  cli::cli_alert_success("Found optimal embedding dimension E to be {E_best}.")
  
  # Find the maximum library size using this E selection
  res <- edm(t, x, E=E_best, full=TRUE, saveManifolds=TRUE, verbosity=verbosity)
  libraryMax <- length(res$Ms[[1]])
  cli::cli_alert_info("The maximum library size we can use is {libraryMax}.")
  
  # Next do causal cross-mapping (CCM) from the cause to the effect
  libraries <- ceiling(seq(10, libraryMax, length.out=25))
  
  res <- edm(t, x, y, E=E_best, library=libraries, verbosity=verbosity)
  
  return()
}