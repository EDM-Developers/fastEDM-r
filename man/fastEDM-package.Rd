\name{fastEDM-package}
\alias{fastEDM-package}
\alias{fastEDM}
\docType{package}
\title{
\packageTitle{fastEDM}
}
\description{
\packageDescription{fastEDM}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{fastEDM}
\packageIndices{fastEDM}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\packageAuthor{fastEDM}

Maintainer: \packageMaintainer{fastEDM}
}
\references{
Jinjing Li, Michael J. Zyphur, George Sugihara, Patrick J. Laub (2021), Beyond Linearity, Stability, and Equilibrium: The edm Package for Empirical Dynamic Modeling and Convergent Cross Mapping in Stata, Stata Journal, 21(1), pp. 220-258
}

\keyword{ package }
\seealso{
\url{https://jinjingli.github.io/edm/}
}
\examples{
\dontrun{
t <- c(1, 2, 3, 4)
x <- c(11, 12, 13, 14)
edm(t, x)
}
}
