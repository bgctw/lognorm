\name{getParmsLognormForLowerAndUpper}
\alias{getParmsLognormForLowerAndUpper}
\title{getParmsLognormForLowerAndUpper}
\description{Calculates mu and sigma of lognormal from lower and upper quantile.}
\usage{getParmsLognormForLowerAndUpper(lower, upper, 
    sigmaFac = qnorm(0.99), isTransScale = FALSE)}
\arguments{
  \item{lower}{value at the lower quantile, i.e. practical minimum}
  \item{upper}{value at the upper quantile, i.e. practical maximum}
  \item{sigmaFac}{sigmaFac = 2 is 95\%
sigmaFac = 2.6 is 99\% interval}
  \item{isTransScale}{if true lower and upper are already on log scale}
}

\value{named numeric vector: mu and sigma parameter of the lognormal
distribution.}

\author{Thomas Wutzler}




\examples{
# sample in normal space
mu <- 5
sigma <- 2
rrNorm <- rnorm(1000, mean = mu, sd = sigma)
# transform to orignal scale
rrOrig <- exp(rrNorm)
# and re-estimate parameters from original scale
res <- getParmsLognormForMedianAndUpper(
  median(rrOrig), quantile(rrOrig, probs = 0.95), sigmaFac = qnorm(0.95))
expected <- c(mu = mu,sigma = sigma)
all.equal(res[1,], expected, tolerance = .1, scale = 1)
}
