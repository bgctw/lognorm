\name{varEffective}
\alias{varEffective}
\title{varEffective}
\description{estimate the variance of a correlated time series}
\usage{varEffective(res, nEff = computeEffectiveNumObs(res), 
    ...)}
\arguments{
  \item{res}{numeric of autocorrelated numbers, usually observation - model residuals}
  \item{nEff}{effective number of observations}
  \item{\dots}{further arguments to \code{\link{var}}}
}
\details{The BLUE is not anymore the usual variance, but a modified
variance as given in Zieba 2011}
\value{The estimated variance of the sample}

\author{Thomas Wutzler}





