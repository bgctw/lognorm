\name{computeEffectiveNumObs}
\alias{computeEffectiveNumObs}
\title{computeEffectiveNumObs}
\description{compute the effective number of observations taking into account autocorrelation}
\usage{computeEffectiveNumObs(res, effAcf = computeEffectiveAutoCorr(res))}
\arguments{
  \item{res}{numeric of autocorrelated numbers, usually observation - model residuals}
  \item{effAcf}{may provide precomputed for efficiency}
}

\value{integer scalar: effective number of observations}

\author{Thomas Wutzler}





