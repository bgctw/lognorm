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
 \item{ Expected value and variance: \code{\link{getLognormMoments}} }
 \item{ Mode: \code{\link{getLognormMode}} }
 \item{ Median: \code{\link{getLognormMedian}} }
 }
    
Estimating parameters 
 \itemize{
 \item{from sample: 
   \code{\link{estimateParmsLognormFromSample}} }
 \item{corresponding to mean and variance at original scale: 
   \code{\link{getParmsLognormForMoments}} }
 \item{corresponding to mean and multiplicative standard deviation at original scale: 
   \code{\link{getParmsLognormForExpval}} }
 }
 
Approximate the sum of correlated lognormals 
 \itemize{
 \item{According to Lo 2013: \code{\link{estimateSumLognormal}} }
 }

Have a look at the package vignettes.
 
}%details

\references{
Limpert E, Stahel W & Abbt M (2001) Log-normal Distributions across the Sciences: 
Keys and Clues.BioScience, Oxford University Press (OUP), 51 , 341 
10.1641/0006-3568(2001)051[0341:lndats]2.0.co;2

Lo C (2013) WKB approximation for the sum of two correlated lognormal 
random variables.
Applied Mathematical Sciences, Hikari, Ltd., 7 , 6355-6367 
10.12988/ams.2013.39511
}

\author{Thomas Wutzler}

\keyword{ package }


