\name{computeEffectiveNumObs}
\alias{computeEffectiveNumObs}
\title{computeEffectiveNumObs}
\description{compute the effective number of observations taking into account autocorrelation}
\usage{computeEffectiveNumObs(res, effAcf = computeEffectiveAutoCorr(res), 
    na.rm = FALSE)}
\arguments{
  \item{res}{numeric of autocorrelated numbers, usually observation - model residuals}
  \item{effAcf}{may provide precomputed for efficiency}
  \item{na.rm}{a logical value indicating whether NA values should be
stripped before the computation proceeds. }
}

\value{integer scalar: effective number of observations}

\author{Thomas Wutzler}





