\name{computeEffectiveAutoCorr}
\alias{computeEffectiveAutoCorr}
\title{computeEffectiveAutoCorr}
\description{return the vector of effective components of the autocorrelation}
\usage{computeEffectiveAutoCorr(res)}
\arguments{
  \item{res}{numeric of autocorrelated numbers, usually observation - model residuals}
}
\details{Returns all components before first negative autocorrelation}
\value{numeric vector: stongest compponents of the autocorrelation function}
\references{Zieba 2011 Standard Deviation of the Mean of Autocorrelated
Observations Estimated with the Use of the Autocorrelation Function Estimated
From the Data}
\author{Thomas Wutzler}





