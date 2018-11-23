\name{computeEffectiveNumObs}
\alias{computeEffectiveNumObs}
\title{computeEffectiveNumObs}
\description{Compute the effective number of observations taking into account autocorrelation}
\usage{computeEffectiveNumObs(res, effAcf = computeEffectiveAutoCorr(res), 
    na.rm = FALSE)}
\arguments{
  \item{res}{numeric of autocorrelated numbers, usually observation - model residuals}
  \item{effAcf}{autocorrelation coefficients.
The first entry is fixed at 1 for zero distance.
May provide precomputed for efficiency or computed from a larger sample.}
  \item{na.rm}{a logical value indicating whether NA values should be
stripped before the computation proceeds. }
}
\details{Handling of NA values: NAs at the beginning or end and are 
just trimmed before computation and pose no problem. 
However with NAs aside from edges, the return value is biased low,
because correlation terms are subtracted for those positions.

Because of NA correlation terms, the computed effective number of
observations can be smaller than 1. In this case 1 is returned.}
\value{integer scalar: effective number of observations}
\references{\code{Zieba & Ramza (2011) 
Standard Deviation of the Mean of Autocorrelated 
Observations Estimated with the Use of the Autocorrelation Function 
Estimated From the Data. 
Metrology and Measurement Systems, 
Walter de Gruyter GmbH, 18 10.2478/v10178-011-0052-x}

\code{Bayley & Hammersley (1946) 
The "effective" number of independent observations in an autocorrelated 
time series. 
Supplement to the Journal of the Royal Statistical Society, JSTOR,8,184-197}}
\author{Thomas Wutzler}




\examples{
# generate autocorrelated time series
res <- stats::filter(rnorm(1000), filter = rep(1,5), circular = TRUE)
res[100:120] <- NA
# plot the series of autocorrelated random variables
plot(res)
# plot their empirical autocorrelation function
acf(res, na.action = na.pass)
#effAcf <- computeEffectiveAutoCorr(res)
# the effective number of parameters is less than number of 1000 samples
(nEff <- computeEffectiveNumObs(res, na.rm = TRUE))
}
