\name{lognorm-package}
\alias{lognorm-package}
\alias{lognorm}
\title{Utilities for the lognormal distribution in R}
\description{
Utilities for the lognormal distribution in R
\itemize{
\item Compute moments.
\item Estimate autocorrelation.
\item Approximate the sum of correlated lognormals.
}
}%description


\details{

Moments and mode 
 \itemize{
 \item Expected value and variance: \code{\link{getLognormMoments}} 
 \item Mode: \code{\link{getLognormMode}} 
 \item Median: \code{\link{getLognormMedian}} 
 }
    
Estimating parameters 
\itemize{
 \item from sample: 
   \code{\link{estimateParmsLognormFromSample}} 
 \item from mean and variance at original scale: 
   \code{\link{getParmsLognormForMoments}} 
 \item from mean and multiplicative standard deviation at original 
   scale: \code{\link{getParmsLognormForExpval}} 
 \item from expected value and upper quantile: 
   \code{\link{getParmsLognormForMeanAndUpper}} 
 \item from median and upper quantile: 
   \code{\link{getParmsLognormForMedianAndUpper}} 
 \item from mode and upper quantile: 
   \code{\link{getParmsLognormForModeAndUpper}} 
 \item from lower and upper quantile: 
   \code{\link{getParmsLognormForLowerAndUpper}} 
}
 
Approximate the sum of correlated lognormals 
\itemize{
 \item According to Lo 2013: \code{\link{estimateSumLognormal}} 
}
 
Utilities for correlated data. These functions maybe moved to
a separate package in future.
\itemize{
 \item Estimate standard error of the mean: \code{\link{seCor}} 
 \item Compute the effective number of observations taking into 
   account autocorrelation: \code{\link{computeEffectiveNumObs}} 
 \item Return the vector of effective components of the autocorrelation: 
   \code{\link{computeEffectiveAutoCorr}} 
 \item Estimate the variance of a correlated time series: 
   \code{\link{varEffective}} 
}

Also have a look at the \href{../doc/index.html}{package vignettes}.
 
}%details

\references{
\code{Limpert E, Stahel W & Abbt M (2001) Log-normal Distributions across the Sciences: 
Keys and Clues.BioScience, Oxford University Press (OUP), 51 , 341 
10.1641/0006-3568(2001)051[0341:lndats]2.0.co;2}

\code{Lo C (2013) WKB approximation for the sum of two correlated lognormal 
random variables.
Applied Mathematical Sciences, Hikari, Ltd., 7 , 6355-6367 
10.12988/ams.2013.39511}
}

\author{Thomas Wutzler}

\keyword{ package }


