\name{computeEffectiveNumObs}
\alias{computeEffectiveNumObs}
\title{computeEffectiveNumObs}
\description{compute the effective number of observations taking into account autocorrelation}
\usage{computeEffectiveNumObs(res, effAcf = computeEffectiveAutoCorr(res), 
    na.rm = FALSE)}
\arguments{
  \item{res}{numeric of autocorrelated numbers, usually observation - model residuals}
  \item{effAcf}{may provide precomputed for efficiency
the first entry is fixed at 1 for zero distance}
  \item{na.rm}{a logical value indicating whether NA values should be
stripped before the computation proceeds. }
}
\details{Handling of NA values: NAs at the beginning or end and are 
just trimmed before computation and pose no problem. 
However with NAs aside from edges, the return value isbiased low,
because correlatations terms are subtracted for those positions.

Because of NA correlation terms, the computed effective number of
observations can be smaller than 1. In this case 1 is returned.}
\value{integer scalar: effective number of observations}
\references{Zieba & Ramza (2011) 
Standard Deviation of the Mean of Autocorrelated 
Observations Estimated with the Use of the Autocorrelation Function 
Estimated From the Data. 
Metrology and Measurement Systems, 
Walter de Gruyter GmbH, 18 10.2478/v10178-011-0052-x

Bayley & Hammersley (1946) 
The "effective" number of independent observations in an autocorrelated 
time series. 
Supplement to the Journal of the Royal Statistical Society, JSTOR,8,184-197}
\author{Thomas Wutzler}





