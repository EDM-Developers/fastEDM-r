
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastEDM <img src="man/figures/logo.png" align="right" height="200" alt="logo" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/fastEDM)](https://CRAN.R-project.org/package=fastEDM)
<!-- badges: end -->

The `fastEDM` R package implements a series of *Empirical Dynamic
Modeling* tools that can be used for *causal analysis of time series*
data.

Key features of the package:

-   powered by a fast multi-threaded *C++ backend*,
-   able to process panel data, a.k.a. *multispatial EDM*,
-   able to handle *missing data* using new `dt` algorithms or by
    dropping points.

## Installation

You can install the development version of fastEDM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("EDM-Developers/fastEDM")
```

## Example: Chicago crime levels and temperature

This example, looking at the causal links between Chicago’s temperature
and crime rates, is described in full in our
[paper](https://jinjingli.github.io/edm/edm-wp.pdf):

``` r
library(fastEDM)
library(readr)

data <- url("https://raw.githubusercontent.com/EDM-Developers/EDM/master/test/chicago.csv")

chicago <- read_csv(data, col_types = cols(crime = col_double()))
chicago <- head(chicago, 500) # Just to speed up the example

easy_edm("crime", "temp", data=chicago)
#> v Found optimal embedding dimension E to be 7.
#> i The maximum library size we can use is 3451.
#> [1] 0
```

## Stata Package

This package is an R port of our [EDM Stata
package](https://edm-developers.github.io/EDM/). As both packages share
the same underlying C++ code, their behaviour will be identical. If you
plan to adjust some of the various low-level EDM parameters, check out
the documentation of the Stata package for more details on their options
and behaviours.
