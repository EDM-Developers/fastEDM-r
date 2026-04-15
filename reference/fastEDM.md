# fastEDM

The fastEDM package implements a series of tools that can be used for
empirical dynamic modeling. The core algorithm is written in C++ to
achieve a reasonable execution speed.

## References

Jinjing Li, Michael J. Zyphur, George Sugihara, Patrick J. Laub (2021),
Beyond Linearity, Stability, and Equilibrium: The edm Package for
Empirical Dynamic Modeling and Convergent Cross Mapping in Stata, Stata
Journal, 21(1), pp. 220-258

## See also

<https://edm-developers.github.io/fastEDM-r/>,
<https://edm-developers.github.io/fastEDM-python/>, and
<https://edm-developers.github.io/edm-stata/>

## Author

Patrick Laub \<patrick.laub@gmail.com\>

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
