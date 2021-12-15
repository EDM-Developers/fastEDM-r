
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastEDM <img src="man/figures/logo.png" align="right" height="200" alt="logo" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/fastEDM)](https://CRAN.R-project.org/package=fastEDM)
<!-- badges: end -->

The fastEDM package implements a series of tools that can be used for
empirical dynamic modeling. The core algorithm is written in C++ to
achieve a reasonable execution speed.

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
library(ggplot2)
library(readr)

data <- url("https://raw.githubusercontent.com/EDM-Developers/EDM/master/test/chicago.csv")
chicago <- read_csv(data, col_types = cols(crime = col_double()))

libs <- c(seq(10, 200, 5), seq(210, 1000, 10), seq(1020, 2000, 20),
          seq(2050, 4350, 50), 4365)

res1 <- edm(chicago["t"], chicago["temp"], chicago["crime"],
                       E=7, library=libs, numReps=4, verbosity=0, numThreads=4)

res2 <- edm(chicago["t"], chicago["crime"], chicago["temp"],
                       E=7, library=libs, numReps=4, verbosity=0, numThreads=4)

averaged1 <- aggregate(res1$summary[, c("mae", "rho")], list(res1$summary$library), mean)
averaged2 <- aggregate(res2$summary[, c("mae", "rho")], list(res2$summary$library), mean)
colnames(averaged1)[[1]] <- "library"
colnames(averaged2)[[1]] <- "library"

p <- ggplot() + 
  geom_line(data = averaged1, aes(x = library, y = rho), color = 2) +
  geom_line(data = averaged2, aes(x = library, y = rho), color = 3) +
  xlab("Library") +
  ylab("Rho")

p <- p + 
  geom_point(data = res1$summary, aes(x = library, y = rho), alpha = 0.05, color = 2) +
  geom_point(data = res2$summary, aes(x = library, y = rho), alpha = 0.05, color = 3) +
  geom_smooth(method = "loess", alpha = 0.9)

print(p)
```

<img src="man/figures/README-example-1.png" width="100%" />
