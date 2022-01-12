
easy_edm <- function(cause, effect, time=NULL, data=NULL,
                     direction="oneway") {
  
  # If doing a rigorous check, begin by seeing if the cause & effect
  # variable appear to be non-linear dynamical system outputs, or just
  # random noise.
  #TODO
  
  cli::cli_alert_info("length of vecs are {length(cause)} and {length(effect)}")
  
  # First find out the embedding dimension of the causal variable
  if (is.null(data)) {
    cli::cli_alert_info("Using supplied time series vectors.")
    x <- cause
    y <- effect
    t <- if (is.null(time)) seq(length(cause)) else time
  } else {
    cli::cli_alert_info("Pulling the time series from the supplied dataframe.")
    x <- data[cause]
    y <- data[effect]
    t <- data[time]
  }
  
  res <- edm(t, x, E=seq(2, 15))
  E_best <- which.max(res$summary$rho)
  cli::cli_alert_success("Found optimal embedding dimension E to be {E_best}.")
  
  # Find the maximum library size using this E selection
  res <- edm(t, x, E=E_best, full=TRUE, saveManifolds=TRUE)
  libraryMax <- length(res$Ms[[1]])
  cli::cli_alert_info("The maximum library size we can use is {libraryMax}.")
  
  # Next do causal cross-mapping (CCM) from the cause to the effect
  libraries <- ceiling(seq(10, libraryMax, length.out=25))
  
  res <- edm(t, x, y, E=E_best, library=libraries)
  
  return(res)
}