#' edm
#'
#' This function provides access to the empirical dynamical modelling (EDM)
#' routines (which are implemented efficiently in C++). The interface is quite
#' low-level, so beginners are recommended to instead use `easy_edm`.
#'
#' @param t The time variable
#'
#' @param x The first variable in the causal analysis
#'
#' @param y The second variable in the causal analysis
#'
#' @param panel If the data is panel data, then this variable specifies the
#' panel ID of each observation.
#'
#' @param E This option specifies the number of dimensions \eqn{E} used for the
#' main variable in the manifold reconstruction. If a list of numbers is provided, the command will
#' compute results for all numbers specified. The xmap subcommand only supports a single integer as the
#' option whereas the explore subcommand supports the option as a numlist. The default value for \eqn{E} is
#' 2, but in theory \eqn{E} can range from 2 to almost half of the total sample size. The actual \eqn{E} used in
#' the estimation may be different if additional variables are incorporated. A error message is
#' provided if the specified value is out of range. Missing data will limit the maximum \eqn{E} under the
#' default deletion method.
#'
#' @param tau The tau (or \eqn{\tau}) option allows researchers to specify the \sQuote{time delay},
#' which essentially sorts the data by the multiple \eqn{\tau}. This is done by specifying lagged embeddings
#' that take the form: \eqn{t,t-\tau,…,t-(E-1)\tau}, where the default is \code{tau=1} (i.e., typical lags). However, if
#' \code{tau=2} is set then every-other \eqn{t} is used to reconstruct the attractor and make predictions—this does
#' not halve the observed sample size because both odd and even \eqn{t} would be used to construct the set of
#' embedding vectors for analysis. This option is helpful when data are oversampled (i.e., spaced too
#' closely in time) and therefore very little new information about a dynamic system is added at each
#' occasion. However, the \code{tau} setting is also useful if different dynamics occur at different times
#' scales, and can be chosen to reflect a researcher’s theory-driven interest in a specific time-scale
#' (e.g., daily instead of hourly). Researchers can evaluate whether \eqn{\tau > 1} is required by checking for
#' large autocorrelations in the observed data. Of course, such
#' a linear measure of association may not work well in nonlinear systems and thus researchers can also
#' check performance by examining \eqn{\rho} and MAE at different values of \eqn{\tau}.
#'
#' @param theta Theta (or \eqn{\theta}) is the distance weighting parameter for the
#' local neighbours in the manifold. It is used to detect the nonlinearity of the system in the explore
#' subcommand for S-mapping. Of course, as noted above, for simplex projection and CCM a weight of
#' \code{theta = 1} is applied to neighbours based on their distance, which is reflected in the fact that the
#' default value of \eqn{\theta} is 1. However, this can be altered even for simplex projection or CCM (two cases
#' that we do not cover here). Particularly, values for S-mapping to test for improved predictions as
#' they become more local may include
#'
#' \code{theta = c(0, .00001, .0001, .001, .005, .01, .05, .1, .5, 1, 1.5, 2, 3, 4, 6, 8, 10)}.
#'
#' @param library This option specifies the total library size \eqn{L} used for
#' the manifold reconstruction. Varying the library size is used to estimate the convergence property
#' of the cross-mapping, with a minimum value \eqn{L_{min} = E + 2} and the maximum equal to the total number of
#' observations minus sufficient lags (e.g., in the time-series case without missing data this is \eqn{L_{max}
#' = T + 1 - E)}. An error message is given if the \eqn{L} value is beyond the allowed range. To assess the
#' rate of convergence (i.e., the rate at which \eqn{\rho} increases as \eqn{L} grows), the full range of library
#' sizes at small values of \eqn{L} can be used, such as if \eqn{E = 2} and \eqn{T = 100}, with the setting then perhaps
#' being \code{library = c(seq(4, 25), seq(30, 50, 5), seq(54, 99, 15))}.
# This option is only available with the `xmap` subcommand.
#'
#' @param k This option specifies the number of neighbours used for prediction. When
#' set to 1, only the nearest neighbour is used, but as \eqn{k} increases the next-closest nearest neighbours
#' are included for making predictions. In the case that \eqn{k} is set 0, the number of neighbours used is
#' calculated automatically (typically as \eqn{k = E + 1} to form a simplex around a target), which is the
#' default value. When \eqn{k = \infty} (e.g., \code{k=Inf}), all possible points in the prediction set are used (i.e.,
#' all points in the library are used to reconstruct the manifold and predict target vectors). This
#' latter setting is useful and typically recommended for S-mapping because it allows all points in the
#' library to be used for predictions with the weightings in theta. However, with large datasets this
#' may be computationally burdensome and therefore \code{k=100} or perhaps \code{k=500} may be preferred if \eqn{T} or \eqn{NT}
#' is large.
#'
#' @param algorithm This option specifies the algorithm used for prediction. If not
#' specified, simplex projection (a locally weighted average) is used. Valid options include simplex
#' and smap, the latter of which is a sequential locally weighted global linear mapping (or S-map as
#' noted previously). In the case of the xmap subcommand where two variables predict each other, the
#' \code{algorithm="smap"} invokes something analogous to a distributed lag model with \eqn{E + 1} predictors
#' (including a constant term \eqn{c}) and, thus, \eqn{E + 1} locally-weighted coefficients for each predicted
#' observation/target vector—because each predicted observation has its own type of regression done
#' with \eqn{k} neighbours as rows and \eqn{E + 1} coefficients as columns. As noted below, in this case special
#' options are available to save these coefficients for post-processing but, again, it is not actually
#' a regression model and instead should be seen as a manifold.
#'
#' @param p This option adjusts the default number of observations ahead
#' which we predict.
# By default, the explore mode predict \eqn{\tau} observations ahead and the xmap mode uses p(0).
#' This parameter can be negative.
#'
#' @param crossfold This option asks the program to run a cross-fold validation of the
#' predicted variables. crossfold(5) indicates a 5-fold cross validation. Note that this cannot be used
#' together with replicate.
# This option is only available with the `explore` subcommand.
#'
#' @param full When this option is specified, the explore command will use all possible
#' observations in the manifold construction instead of the default 50/50 split. This is effectively
#' the same as leave-one-out cross-validation as the observation itself is not used for the prediction.
#'
#' @param shuffle When splitting the observations into library and prediction sets, by default
#' the oldest observations go into the library set and the newest observations to the prediction set.
#' Though if the randomize option is specified, the data is allocated into the two sets in a random
#' fashion. If the replicate option is specified, then this randomization is enabled automatically.
#'
#' @param copredict This option specifies the variable used for coprediction.
#' A second prediction is run for each configuration of \eqn{E}, library, etc., using the same library set
#' but with a prediction set built from the lagged embedding of this variable.
#'
#' @param savePredictions This option allows you to save the edm predictions
#' which could be useful for plotting and diagnosis.
#'
#' @param saveCoPredictions This option allows you to save the copredictions.
#' You must specify the \code{copredict} option for this to work.
#'
#' @param saveManifolds This option allows you to save the library and prediction manifolds.
#'
#' @param saveSMAPCoeffs  This option allows S-map coefficients to be stored.
# in variables with
# a specified prefix. For example, specifying "edm xmap x y, algorithm(smap) savesmap(beta) k(-1)" will
# create a set of new variables such as beta1_b0_rep1. The string prefix (e.g., \sQuote{beta}) must not be
# shared with any variables in the dataset, and the option is only valid if the algorithm(smap) is
# specified. In terms of the saved variables such as beta1_b0_rep1, the first number immediately after
# the prefix \sQuote{beta} is 1 or 2 and indicates which of the two listed variables is treated as the
# dependent variable in the cross-mapping (i.e., the direction of the mapping). For the "edm xmap x y"
# case, variables starting with beta1_ contain coefficients derived from the manifold \eqn{M_X} created
# using the lags of the first variable 'x' to predict \eqn{Y}, or \eqn{Y|M_X}. This set of variables therefore
# store the coefficients related to 'x' as an outcome rather than a predictor in CCM. Keep in mind
# that any \eqn{Y \to X} effect associated with the beta1_ prefix is shown as \eqn{Y|M_X}, because the outcome is used
# to cross-map the predictor, and thus the reported coefficients will be scaled in the opposite
# direction of a typical regression (because in CCM the outcome variable predicts the cause). To get
# more familiar regression coefficients (which will be locally weighted), variables starting with
# beta2_ store the coefficients estimated in the other direction, where the second listed variable \sQuote{y}
# is used for the manifold reconstruction \eqn{M_Y} for the mapping \eqn{X|M_Y} in the "edm xmap x y" case, testing
# the opposite \eqn{X \to Y} effect in CCM, but with reported S-map coefficients that map to a \eqn{Y \to X} regression.
# We appreciate that this may be unintuitive, but because CCM causation is tested by predicting the
# causal variable with the outcome, to get more familiar regression coefficients requires reversing
# CCM’s causal direction to a more typical predictor -> outcome regression logic. This can be clarified
# by reverting to the conditional notation such as \eqn{X|M_Y}, which in CCM implies a left-to-right \eqn{X \to Y}
# effect, but for the S-map coefficients will be scaled as a locally-weighted regression in the
# opposite direction \eqn{Y \to X}. Moving on, following the 1 and 2 is the letter b and a number. The numerical
# labeling scheme generally follows the order of the lag for the main variable and then the order of
# the extra variables introduced in the case of multivariate embedding. b0 is a special case which
# records the coefficient of the constant term in the regression. The final term rep1 indicates the
# coefficients are from the first round of replication (if the replicate() option is not used then
# there is only one). Finally, the coefficients are saved to match the observation \eqn{t} in the dataset
# that is being predicted, which allows plotting each of the \eqn{E} estimated coefficients against time
# and/or the values of the variable being predicted. The variables are also automatically labelled for
# clarity. This option is only available with the `xmap` subcommand.
#'
#' @param extras This option allows incorporating additional variables into the
#' embedding (multivariate embedding), e.g. \code{extras=c(y, z)}.
# Time series lists are unabbreviated here,
# e.g. extra(L(1/3).z) will be equivalent to extra(L1.z L2.z L3.z). Normally, lagged versions of the
# extra variables are not included in the embedding, however the syntax extra(z(e)) includes e lags
# of z in the embedding.
#'
#' @param allowMissing This option allows observations with missing values to be used in the
#' manifold. Vectors with at least one non-missing values will be used in the manifold construction.
#' Distance computations are adapted to allow missing values when this option is specified.
#'
#' @param missingDistance This option allows users to specify the assumed distance
#' between missing values and any values (including missing) when estimating the Euclidean distance of
#' the vector. This enables computations with missing values. The option implies \code{allowmissing}. By
#' default, the distance is set to the expected distance of two random draws in a normal distribution,
#' which equals to \eqn{2/\sqrt{\pi} \times} standard deviation of the mapping variable.
#'
#' @param dt This option allows automatic inclusion of the timestamp differencing in the
#' embedding. There will be \eqn{E} dt variables included for an embedding with \eqn{E} dimensions. By default,
#' the weights used for these additional variables equal to the standard deviation of the main
#' mapping variable divided by the standard deviation of the time difference. This can be overridden by
#' the \code{dtWeight} option. The \code{dt} option will be ignored when running with data with no sampling
#' variation in the time lags. The first dt variable embeds the time of the between the most recent
#' observation and the time of the corresponding target/predictand.
#'
#' @param reldt This option, to be read as \sQuote{relative dt}, is like the \code{dt} option above in
#' that it includes \eqn{E} extra variables for an embedding with E dimensions. However the timestamp
#' differences added are not the time between the corresponding observations, but the time of the
#' target/predictand minus the time of the lagged observations.
#'
#' @param dtWeight This option specifies the weight used for the timestamp differencing
#' variable.
#'
#' @param numReps The number of random replications (i.e. random splits to library and prediction sets) to run.
# The explore subcommand uses a random 50/50 split for simplex
# projection and S-maps, whereas the xmap subcommand selects the observations randomly for library
# construction if the size of the library \eqn{L} is smaller than the size of all available observations. In
# these cases, results may be different in each run because the embedding vectors (i.e., the
# E-dimensional points) used to reconstruct a manifold are chosen at random.
# The replicate option
# takes advantages of this to allow repeating the randomization process and calculating results each time.
#' This is akin to a nonparametric bootstrap without replacement, and is commonly used for
#' inference using confidence intervals in EDM (Tsonis et al., 2015; van Nes et al., 2015; Ye et al.,
#' 2015b).
# When replicate is specified, such as replicate(50), mean values and the standard deviations
# of the results are reported across the 50 runs by default. As we note below, it is possible to save
# all estimates for post-processing using typical Stata commands such as svmat, allowing the graphing
# of results or finding percentile-based with the pctile command.
#'
#' @param panelWeight This specifies a penalty that is added to the distances between points in the
#' manifold which correspond to observations from different panels. By default \code{panelWeight} is 0,
#' so the data from all panels is mixed together and treatly equally. If \code{panelWeight=Inf} is set
#' then the weight is treated as \eqn{\infty} so neighbours will never
#' be selected which cross the boundaries between panels. Setting \code{panelWeight=Inf} with
#' \code{k=Inf} means we may use a different number of neighbors for different predictions (i.e.
#' if the panels are unbalanced).
#'
#' @param panelWeights A generalisation of \code{panelWeight}. Instead of giving
#' a constant penalty for differing panels, \code{panelWeights} lets you supply a
#' matrix so \code{panelWeights[i, j]} will be added to distances between points
#' in the $i$-th panel and the $j$-th panel.
#'
#' @param verbosity The level of detail in the output.
#'
#' @param showProgressBar Whether or not to print out a progress bar during the computations.
#'
#' @param numThreads The number of threads to use for the prediction task.
#'
#' @param lowMemory The lowMemory option tries to save as much space as possible
#' by more efficiently using memory, though for small datasets this will likely
#' slow down the computations by a small but noticeable amount.
#'
#' @param predictWithPast Force all predictions to only use contemporaneous
#' data. Normally EDM is happy to cheat by pulling segments from the future
#' of the time series to make a prediction.
#'
#' @return A list
#' @export
#'
#' @examples
#' t <- c(1, 2, 3, 4, 5, 6, 7, 8)
#' x <- c(11, 12, 13, 14, 15, 16, 17, 18)
#' res <- edm(t, x)
#' #
edm <- function(t, x, y = c(), panel = c(), E = 2, tau = 1, theta = 1, library = NULL,
                k = 0, algorithm = "simplex", p = NULL, crossfold = 0, full = FALSE,
                shuffle = FALSE, copredict = c(), savePredictions = FALSE,
                saveCoPredictions = FALSE, saveManifolds = FALSE,
                saveSMAPCoeffs = FALSE, extras = NULL, allowMissing = FALSE,
                missingDistance = 0.0, dt = FALSE, reldt = FALSE, dtWeight = 0.0,
                numReps = 1, panelWeight = 0, panelWeights = NULL, verbosity = 1,
                showProgressBar = NA, numThreads = 1, lowMemory = FALSE,
                predictWithPast = FALSE) {
  if (length(t) != length(x)) {
    stop("The time and x variables should be the same length")
  }

  df <- data.frame(t = t, x = x)
  colnames(df) <- c("t", "x")

  if (length(y) > 0) {
    if (length(t) != length(y)) {
      stop("The y variable is the wrong length")
    }
    df["y"] <- y
  }

  if (length(panel) > 0) {
    if (length(t) != length(panel)) {
      stop("The panel id variable is the wrong length")
    }
    df["panel"] <- panel
  }

  if (length(copredict) > 0) {
    if (length(t) != length(copredict)) {
      stop("Coprediction vector is the wrong length")
    }
    df["co_x"] <- copredict
  }

  if (!is.null(extras)) {
    for (extra in extras) {
      if (length(extra) != length(x)) {
        stop("An extra variable is not the right size")
      }
    }
  }

  if (numReps > 1) {
    shuffle <- TRUE
  }

  if (isnothing(showProgressBar)) {
    showProgressBar <- (verbosity > 0)
  }

  explore <- length(y) == 0

  # Re-assert default arguments if NA/NULL/NaN are passed to them
  p <- if (is.null(p)) {
    explore
  } else {
    p
  }
  k <- if (isnothing(k)) {
    0
  } else {
    k
  }
  if (length(E) == 1) {
    E <- if (isnothing(E)) {
      2
    } else {
      E
    }
  } else {
    E <- if (sum(is.na(E)) > 0) {
      E[!is.na(E)]
    } else {
      E
    }
  }

  k <- if (k == Inf) {
    -1
  } else {
    k
  }

  panelWeight <- if (panelWeight == Inf) {
    -1
  } else {
    panelWeight
  }

  res <- run_command(df, E, tau, theta, library,
    k,
    algorithm = algorithm, numReps = numReps,
    p = p, crossfold = crossfold, full = full, shuffle = shuffle,
    saveFinalPredictions = savePredictions,
    saveFinalCoPredictions = saveCoPredictions,
    saveManifolds = saveManifolds,
    saveSMAPCoeffs = saveSMAPCoeffs,
    extras = extras, allowMissing = allowMissing,
    missingDistance = missingDistance,
    dt = dt, reldt = reldt, dtWeight = dtWeight,
    numThreads = numThreads, panelWeight = panelWeight,
    panelWeights = panelWeights, verbosity = verbosity,
    showProgressBar = showProgressBar, lowMemory = lowMemory,
    predictWithPast = predictWithPast
  )

  if (res$rc == 0) {
    if (verbosity > 1) {
      cat("Finished successfully!\n")
    }

    df <- stats::na.omit(res$stats)

    if (verbosity > 1) {
      cat("Number of non-missing stats:", nrow(df), "\n")
    }

    if (nrow(df) > 1) {
      res$summary <- stats::aggregate(cbind(rho, mae) ~ E + library + theta, df, mean)
    } else {
      res$summary <- res$stats
    }


    if (verbosity > 0) {
      cat("Summary of predictions\n")
      print.data.frame(res$summary)

      if (!is.na(res$kMin) && !is.na(res$kMax)) {
        if (res$kMin == res$kMax) {
          cat("Number of neighbours (k) is set to", res$kMin, "\n")
        } else {
          cat(
            "Number of neighbours (k) is set to between",
            res$kMin, "and", res$kMax, "\n"
          )
        }
      }
    }

    if (length(copredict) > 0) {
      df <- stats::na.omit(res$copredStats)

      if (nrow(df) > 1) {
        res$copredSummary <- stats::aggregate(cbind(rho, mae) ~ E + library + theta, df, mean)
      } else {
        res$copredSummary <- res$copredStats
      }

      if (verbosity > 0) {
        cat("Summary of copredictions\n")
        print.data.frame(res$copredSummary)
      }
    }
  } else {
    cat("Error code:", res$rc, "\n")
  }

  return(res)
}

isnothing <- function(x) {
  is.null(x) | is.na(x) | is.nan(x)
}
