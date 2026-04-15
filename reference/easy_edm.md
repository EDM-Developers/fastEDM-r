# easy_edm

This is an automated workflow for performing causal analysis on the
supplied time series using empirical dynamical modelling (EDM)
techniques. It is intended to hide all the common steps of an EDM
analysis, and should work on most datasets.

## Usage

``` r
easy_edm(
  cause,
  effect,
  time = NULL,
  data = NULL,
  direction = "oneway",
  verbosity = 1,
  showProgressBar = NULL,
  normalize = TRUE
)
```

## Arguments

- cause:

  The causal time series (as a string or a vector).

- effect:

  The effect time series (as a string or a vector).

- time:

  For non-regularly sampled time series, the sampling times must be
  supplied here (as a string or a vector).

- data:

  If a dataframe is supplied here, then the cause, effect & time
  arguments must be the column names of the relevant time series as
  strings.

- direction:

  A string specifying whether we are checking a one directional causal
  effect or whether to test the reverse direction at the same time (work
  in progress!).

- verbosity:

  The level of detail in the output.

- showProgressBar:

  Whether or not to print out a progress bar during the computations.

- normalize:

  Whether to normalize the inputs before starting EDM.

## Value

A Boolean indicating that evidence of causation was found.

## Details

It may be the case that your data requires a custom analysis, so this
function can be used as a helpful starting point from which to create a
specialised analysis using the \`edm\` function directly.

Warning: While the \`edm\` functionality is well-tested and ready for
use, this \`easy_edm\` automated analysis is still a work-in-progress.

## Examples

``` r
library(fastEDM)

num_obs <- 500
map <- logistic_map(num_obs)
df <- data.frame(list(x = map$x, y = map$y))

x_causes_y <- easy_edm("x", "y", data=df, verbosity=0)
#> ✔ Strong evidence of CCM causation from x to y found.
y_causes_x <- easy_edm("y", "x", data=df, verbosity=0)
#> ✔ Strong evidence of CCM causation from y to x found.
```
