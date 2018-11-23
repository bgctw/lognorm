\name{varEffective}
\alias{varEffective}
\title{varEffective}
\description{Estimate the variance of a correlated time series}
\usage{varEffective(res, nEff = computeEffectiveNumObs(res, 
    na.rm = na.rm), na.rm = FALSE, ...)}
\arguments{
  \item{res}{numeric of autocorrelated numbers, usually observation - model residuals}
  \item{nEff}{effective
number of observations}
  \item{na.rm}{set to TRUE to remove NA cases before computation}
  \item{\dots}{further arguments to \code{\link{var}}}
}
\details{The BLUE is not anymore the usual variance, but a modified
variance as given in \code{Zieba 2011}}
\value{The estimated variance of the sample}

\author{Thomas Wutzler}




\examples{
# generate autocorrelated time series
res <- stats::filter(rnorm(1000), filter = rep(1,5), circular = TRUE)
res[100:120] <- NA
# if correlations are neglected, the estimate of the variance is biased low
(varNeglectCorr <- var(res, na.rm = TRUE))
(varCorr <- varEffective(res, na.rm = TRUE))
}
