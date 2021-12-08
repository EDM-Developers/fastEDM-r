
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastEDM <img src="man/figures/logo.png" align="right" height="200" alt="logo" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/fastEDM)](https://CRAN.R-project.org/package=fastEDM)
<!-- badges: end -->

The goal of fastEDM is to …

## Installation

You can install the development version of fastEDM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("EDM-Developers/fastEDM")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(fastEDM)

t <- c(1, 2, 3, 4, 5, 6, 7, 8)
x <- c(11, 12, 13, 14, 15, 16, 17, 18)
res <- edm(t, x)
res$summary
#>   E library theta        rho     mae
#> 1 2       3     1 -0.9706895 2.62881
```
