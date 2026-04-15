# edm

This function provides access to the empirical dynamical modelling (EDM)
routines (which are implemented efficiently in C++). The interface is
quite low-level, so beginners are recommended to instead use
\`easy_edm\`.

## Usage

``` r
edm(
  t,
  x,
  y = c(),
  panel = c(),
  E = 2,
  tau = 1,
  theta = 1,
  library = NULL,
  k = 0,
  algorithm = "simplex",
  p = NULL,
  crossfold = 0,
  full = FALSE,
  shuffle = FALSE,
  seed = 0,
  copredict = c(),
  savePredictions = FALSE,
  saveCoPredictions = FALSE,
  saveManifolds = FALSE,
  saveSMAPCoeffs = FALSE,
  extras = NULL,
  allowMissing = FALSE,
  missingDistance = 0,
  dt = FALSE,
  reldt = FALSE,
  dtWeight = 0,
  numReps = 1,
  panelWeight = 0,
  panelWeights = NULL,
  verbosity = 1,
  showProgressBar = NA,
  numThreads = 1,
  lowMemory = FALSE,
  predictWithPast = FALSE
)
```

## Arguments

- t:

  The time variable

- x:

  The first variable in the causal analysis

- y:

  The second variable in the causal analysis

- panel:

  If the data is panel data, then this variable specifies the panel ID
  of each observation.

- E:

  This option specifies the number of dimensions \\E\\ used for the main
  variable in the manifold reconstruction. If a list of numbers is
  provided, the command will compute results for all numbers specified.
  The xmap subcommand only supports a single integer as the option
  whereas the explore subcommand supports the option as a numlist. The
  default value for \\E\\ is 2, but in theory \\E\\ can range from 2 to
  almost half of the total sample size. The actual \\E\\ used in the
  estimation may be different if additional variables are incorporated.
  A error message is provided if the specified value is out of range.
  Missing data will limit the maximum \\E\\ under the default deletion
  method.

- tau:

  The tau (or \\\tau\\) option allows researchers to specify the ‘time
  delay’, which essentially sorts the data by the multiple \\\tau\\.
  This is done by specifying lagged embeddings that take the form:
  \\t,t-\tau,…,t-(E-1)\tau\\, where the default is `tau=1` (i.e.,
  typical lags). However, if `tau=2` is set then every-other \\t\\ is
  used to reconstruct the attractor and make predictions—this does not
  halve the observed sample size because both odd and even \\t\\ would
  be used to construct the set of embedding vectors for analysis. This
  option is helpful when data are oversampled (i.e., spaced too closely
  in time) and therefore very little new information about a dynamic
  system is added at each occasion. However, the `tau` setting is also
  useful if different dynamics occur at different times scales, and can
  be chosen to reflect a researcher’s theory-driven interest in a
  specific time-scale (e.g., daily instead of hourly). Researchers can
  evaluate whether \\\tau \> 1\\ is required by checking for large
  autocorrelations in the observed data. Of course, such a linear
  measure of association may not work well in nonlinear systems and thus
  researchers can also check performance by examining \\\rho\\ and MAE
  at different values of \\\tau\\.

- theta:

  Theta (or \\\theta\\) is the distance weighting parameter for the
  local neighbours in the manifold. It is used to detect the
  nonlinearity of the system in the explore subcommand for S-mapping. Of
  course, as noted above, for simplex projection and CCM a weight of
  `theta = 1` is applied to neighbours based on their distance, which is
  reflected in the fact that the default value of \\\theta\\ is 1.
  However, this can be altered even for simplex projection or CCM (two
  cases that we do not cover here). Particularly, values for S-mapping
  to test for improved predictions as they become more local may include

  `theta = c(0, .00001, .0001, .001, .005, .01, .05, .1, .5, 1, 1.5, 2, 3, 4, 6, 8, 10)`.

- library:

  This option specifies the total library size \\L\\ used for the
  manifold reconstruction. Varying the library size is used to estimate
  the convergence property of the cross-mapping, with a minimum value
  \\L\_{min} = E + 2\\ and the maximum equal to the total number of
  observations minus sufficient lags (e.g., in the time-series case
  without missing data this is \\L\_{max} = T + 1 - E)\\. An error
  message is given if the \\L\\ value is beyond the allowed range. To
  assess the rate of convergence (i.e., the rate at which \\\rho\\
  increases as \\L\\ grows), the full range of library sizes at small
  values of \\L\\ can be used, such as if \\E = 2\\ and \\T = 100\\,
  with the setting then perhaps being
  `library = c(seq(4, 25), seq(30, 50, 5), seq(54, 99, 15))`.

- k:

  This option specifies the number of neighbours used for prediction.
  When set to 1, only the nearest neighbour is used, but as \\k\\
  increases the next-closest nearest neighbours are included for making
  predictions. In the case that \\k\\ is set 0, the number of neighbours
  used is calculated automatically (typically as \\k = E + 1\\ to form a
  simplex around a target), which is the default value. When \\k =
  \infty\\ (e.g., `k=Inf`), all possible points in the prediction set
  are used (i.e., all points in the library are used to reconstruct the
  manifold and predict target vectors). This latter setting is useful
  and typically recommended for S-mapping because it allows all points
  in the library to be used for predictions with the weightings in
  theta. However, with large datasets this may be computationally
  burdensome and therefore `k=100` or perhaps `k=500` may be preferred
  if \\T\\ or \\NT\\ is large.

- algorithm:

  This option specifies the algorithm used for prediction. If not
  specified, simplex projection (a locally weighted average) is used.
  Valid options include simplex and smap, the latter of which is a
  sequential locally weighted global linear mapping (or S-map as noted
  previously). In the case of the xmap subcommand where two variables
  predict each other, the `algorithm="smap"` invokes something analogous
  to a distributed lag model with \\E + 1\\ predictors (including a
  constant term \\c\\) and, thus, \\E + 1\\ locally-weighted
  coefficients for each predicted observation/target vector—because each
  predicted observation has its own type of regression done with \\k\\
  neighbours as rows and \\E + 1\\ coefficients as columns. As noted
  below, in this case special options are available to save these
  coefficients for post-processing but, again, it is not actually a
  regression model and instead should be seen as a manifold.

- p:

  This option adjusts the default number of observations ahead which we
  predict. This parameter can be negative.

- crossfold:

  This option asks the program to run a cross-fold validation of the
  predicted variables. crossfold(5) indicates a 5-fold cross validation.
  Note that this cannot be used together with replicate.

- full:

  When this option is specified, the explore command will use all
  possible observations in the manifold construction instead of the
  default 50/50 split. This is effectively the same as leave-one-out
  cross-validation as the observation itself is not used for the
  prediction.

- shuffle:

  When splitting the observations into library and prediction sets, by
  default the oldest observations go into the library set and the newest
  observations to the prediction set. Though if the randomize option is
  specified, the data is allocated into the two sets in a random
  fashion. If the replicate option is specified, then this randomization
  is enabled automatically.

- seed:

  To allow for reproducibility of the internal random number generator,
  set this seed to a given integer value.

- copredict:

  This option specifies the variable used for coprediction. A second
  prediction is run for each configuration of \\E\\, library, etc.,
  using the same library set but with a prediction set built from the
  lagged embedding of this variable.

- savePredictions:

  This option allows you to save the edm predictions which could be
  useful for plotting and diagnosis.

- saveCoPredictions:

  This option allows you to save the copredictions. You must specify the
  `copredict` option for this to work.

- saveManifolds:

  This option allows you to save the library and prediction manifolds.

- saveSMAPCoeffs:

  This option allows S-map coefficients to be stored.

- extras:

  This option allows incorporating additional variables into the
  embedding (multivariate embedding), e.g. `extras=c(y, z)`.

- allowMissing:

  This option allows observations with missing values to be used in the
  manifold. Vectors with at least one non-missing values will be used in
  the manifold construction. Distance computations are adapted to allow
  missing values when this option is specified.

- missingDistance:

  This option allows users to specify the assumed distance between
  missing values and any values (including missing) when estimating the
  Euclidean distance of the vector. This enables computations with
  missing values. The option implies `allowmissing`. By default, the
  distance is set to the expected distance of two random draws in a
  normal distribution, which equals to \\2/\sqrt{\pi} \times\\ standard
  deviation of the mapping variable.

- dt:

  This option allows automatic inclusion of the timestamp differencing
  in the embedding. There will be \\E\\ dt variables included for an
  embedding with \\E\\ dimensions. By default, the weights used for
  these additional variables equal to the standard deviation of the main
  mapping variable divided by the standard deviation of the time
  difference. This can be overridden by the `dtWeight` option. The `dt`
  option will be ignored when running with data with no sampling
  variation in the time lags. The first dt variable embeds the time of
  the between the most recent observation and the time of the
  corresponding target/predictand.

- reldt:

  This option, to be read as ‘relative dt’, is like the `dt` option
  above in that it includes \\E\\ extra variables for an embedding with
  E dimensions. However the timestamp differences added are not the time
  between the corresponding observations, but the time of the
  target/predictand minus the time of the lagged observations.

- dtWeight:

  This option specifies the weight used for the timestamp differencing
  variable.

- numReps:

  The number of random replications (i.e. random splits to library and
  prediction sets) to run. This is akin to a nonparametric bootstrap
  without replacement, and is commonly used for inference using
  confidence intervals in EDM (Tsonis et al., 2015; van Nes et al.,
  2015; Ye et al., 2015b).

- panelWeight:

  This specifies a penalty that is added to the distances between points
  in the manifold which correspond to observations from different
  panels. By default `panelWeight` is 0, so the data from all panels is
  mixed together and treatly equally. If `panelWeight=Inf` is set then
  the weight is treated as \\\infty\\ so neighbours will never be
  selected which cross the boundaries between panels. Setting
  `panelWeight=Inf` with `k=Inf` means we may use a different number of
  neighbors for different predictions (i.e. if the panels are
  unbalanced).

- panelWeights:

  A generalisation of `panelWeight`. Instead of giving a constant
  penalty for differing panels, `panelWeights` lets you supply a matrix
  so `panelWeights[i, j]` will be added to distances between points in
  the \$i\$-th panel and the \$j\$-th panel.

- verbosity:

  The level of detail in the output.

- showProgressBar:

  Whether or not to print out a progress bar during the computations.

- numThreads:

  The number of threads to use for the prediction task.

- lowMemory:

  The lowMemory option tries to save as much space as possible by more
  efficiently using memory, though for small datasets this will likely
  slow down the computations by a small but noticeable amount.

- predictWithPast:

  Force all predictions to only use contemporaneous data. Normally EDM
  is happy to cheat by pulling segments from the future of the time
  series to make a prediction.

## Value

A list

## Examples

``` r
t <- c(1, 2, 3, 4, 5, 6, 7, 8)
x <- c(11, 12, 13, 14, 15, 16, 17, 18)
res <- edm(t, x)
#> Computing: [========================================] 100% (done)                         
#> Summary of predictions
#>   E library theta       rho       mae
#> 1 2       7     1 0.9398467 0.6820797
#> Number of neighbours (k) is set to 3 
#
```
