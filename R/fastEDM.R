
#' fastEDM
#' 
#' The fastEDM package implements a series of tools that can be
#' used for empirical dynamic modeling. The core algorithm is written in C++
#' to achieve a reasonable execution speed. 
#' 
#' @docType package
#' @author Patrick Laub <patrick.laub@gmail.com>
#' @import Rcpp RcppEigen
#' @importFrom Rcpp evalCpp
#' @useDynLib fastEDM
#' @name fastEDM
#' @references Jinjing Li, Michael J. Zyphur, George Sugihara, Patrick J. Laub (2021), Beyond Linearity, Stability, and Equilibrium: The edm Package for Empirical Dynamic Modeling and Convergent Cross Mapping in Stata, Stata Journal, 21(1), pp. 220-258
#' @seealso{
#'   \url{https://jinjingli.github.io/edm/}
#' }
#' @examples 
#'\dontrun{
#'  t <- c(1, 2, 3, 4)
#'  x <- c(11, 12, 13, 14)
#'  edm(t, x)
#'}
NULL