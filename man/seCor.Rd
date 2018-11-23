\name{seCor}
\alias{seCor}
\title{seCor}
\description{Compute the standard error accounting for empirical autocorrelations}
\usage{seCor(x, na.rm = FALSE, effCov = computeEffectiveAutoCorr(x, 
    type = "covariance"))}
\arguments{
  \item{x}{numeric vector}
  \item{na.rm}{logical. Should missing values be removed?}
  \item{effCov}{numeric vector of effective covariance components
first entry is the variance. See \code{\link{computeEffectiveAutoCorr}}}
}
\details{Computation follows 
https://stats.stackexchange.com/questions/274635/calculating-error-of-mean-of-time-series.


The default uses empiricial autocorrelation
estimates from the supplied data up to first negative component.
For short series of \code{x} it is strongly recommended to to
provide \code{effCov} that was estimated on a longer time series.}
\value{numeric scalar of standard error of the mean of x}

\author{Thomas Wutzler}





