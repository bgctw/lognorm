\name{getParmsLognormForMedianAndUpper}
\alias{getParmsLognormForMedianAndUpper}
\title{getParmsLognormForMedianAndUpper}
\description{Calculates mu and sigma of lognormal from median and upper quantile.}
\usage{getParmsLognormForMedianAndUpper(median, 
    quant, sigmaFac = qnorm(0.99))}
\arguments{
  \item{median}{geometric mu (median at the original exponential scale)}
  \item{quant}{value at the upper quantile, i.e. practical maximum}
  \item{sigmaFac}{sigmaFac=2 is 95\% sigmaFac=2.6 is 99\% interval}
}

\value{named numeric vector: mu and sigma parameter of the lognormal distribution.}

\author{Thomas Wutzler}





