library(dplyr)

edm <- function(t, x, y = c(), E=2, tau=1, theta=1, library=NULL, k=0,
                algorithm="simplex", p=NULL, crossfold=0, full=FALSE,
                shuffle=FALSE, copredict = c(), savePredictions=FALSE,
                saveCoPredictions=FALSE, saveSMAPCoeffs=FALSE,
                extras=NULL, allowMissing=FALSE, missingDistance=0.0,
                dt=FALSE, reldt=FALSE, dtWeight=0.0, # saveDT=FALSE,
                numReps=1, verbosity=0, numThreads=1, saveInputs="") {
  
  if (length(y) > 0) {
    df <- data.frame(t = t, x = x, y = y)
    colnames(df) <- c("t", "x", "y")
    explore <- FALSE
  } else {
    df <- data.frame(t = t, x = x)
    colnames(df) <- c("t", "x")
    explore <- TRUE
  }
  
  if (length(copredict) > 0) {
    if (length(copredict) != length(x)) {
      stop("Coprediction vector is not the right size")
    }
    df["co_x"] <- copredict
  }
  
  if (numReps > 1) {
    shuffle = TRUE
  }
  
  p <- if (!is.null(p)) { p } else { explore }
  
  res <- run_command(df, E, tau, theta, library,
                     k, algorithm=algorithm, numReps=numReps,
                     p=p, crossfold=crossfold, full=full, shuffle=shuffle,
                     saveFinalPredictions=savePredictions,
                     saveFinalCoPredictions=saveCoPredictions,
                     saveSMAPCoeffs=saveSMAPCoeffs,
                     extras=extras,  allowMissing=allowMissing,
                     missingDistance=missingDistance,
                     dt=dt, reldt=reldt, dtWeight=dtWeight, 
                     numThreads=numThreads, verbosity=verbosity)
  
  return(res)
}
