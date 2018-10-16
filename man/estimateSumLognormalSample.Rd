\name{estimateSumLognormalSample}
\alias{estimateSumLognormalSample}
\title{estimateSumLognormalSample}
\description{Estimate the parameters of the lognormal approximation to the sum}
\usage{estimateSumLognormalSample(mu, sigma, resLog, 
    effAcf = computeEffectiveAutoCorr(resLog), 
    isGapFilled = logical(0), na.rm = TRUE)}
\arguments{
  \item{mu}{numeric vector of center parameters of terms at log scale}
  \item{sigma}{numeric vector of variance parameter of terms at log scale}
  \item{resLog}{time series of model-residuals at log scale
to estimate correlation}
  \item{effAcf}{effective autocorrelation
coefficients (may provide precomputed for efficiency or if the sample
of \code{resLog} is too small) set to 1 to assume uncorrelated sample}
  \item{isGapFilled}{logical vector whether entry is gap-filled
rather than an original measurement, see details}
  \item{na.rm}{neglect terms with NA values in mu or sigma}
}
\details{If there are no gap-filled values, i.e. \code{all(!isGapFilled)} or
\code{!length(isGapFilled)} (the default), distribution parameters
are estimated using all the samples. Otherwise, the scale parameter
(uncertainty) is first estimated using only the non-gapfilled records.

If there are only gap-filled records, assume uncertainty
to be the largest uncertainty of given gap-filled records.}
\value{numeric vector with components \code{mu}, \code{sigma}, 
and \code{nEff},
the parameters of the lognormal distribution at log scale
(Result of \code{link{estimateSumLognormal}})
and the number of effective observations.}

\author{Thomas Wutzler}





