% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autocorr.R
\name{seCor}
\alias{seCor}
\title{Compute the standard error accounting for empirical autocorrelations}
\usage{
seCor(
  x,
  effCor = if (missing(effCov)) computeEffectiveAutoCorr(x) else effCov/var(x, na.rm =
    TRUE),
  na.rm = FALSE,
  effCov,
  nEff = computeEffectiveNumObs(x, effCor, na.rm = na.rm)
)
}
\arguments{
\item{x}{numeric vector}

\item{effCor}{numeric vector of effective correlation components
first entry at zero lag equals one. See \code{\link{computeEffectiveAutoCorr}}}

\item{na.rm}{logical. Should missing values be removed?}

\item{effCov}{alternative to specifying effCor: numeric vector of 
effective covariance components
first entry is the variance. See \code{\link{computeEffectiveAutoCorr}}}

\item{nEff}{possibility to specify precomputed number of effective 
observations for speedup.}
}
\value{
numeric scalar of standard error of the mean of x
}
\description{
Compute the standard error accounting for empirical autocorrelations
}
\details{
The default uses empirical autocorrelation
 estimates from the supplied data up to first negative component.
 For short series of \code{x} it is strongly recommended to to
 provide \code{effCov} that was estimated on a longer time series.
}
