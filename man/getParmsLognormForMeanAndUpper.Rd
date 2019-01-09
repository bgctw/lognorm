\name{getParmsLognormForMeanAndUpper}
\alias{getParmsLognormForMeanAndUpper}
\title{getParmsLognormForMeanAndUpper}
\description{Calculates mu and sigma of lognormal from median and upper quantile.}
\usage{getParmsLognormForMeanAndUpper(mean, quant, 
    sigmaFac = qnorm(0.99))}
\arguments{
  \item{mean}{expected value at the original scale}
  \item{quant}{value at the upper quantile, i.e. practical maximum}
  \item{sigmaFac}{sigmaFac=2 is 95\% sigmaFac=2.6 is 99\% interval}
}
\details{There are two valid solutions. This routine returns
the one with lower sigma, i.e. the not so strongly skewed solution.}
\value{numeric matrix: columns mu and sigma parameter of the
lognormal distribution.}

\author{Thomas Wutzler}





