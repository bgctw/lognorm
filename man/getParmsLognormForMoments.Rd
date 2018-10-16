\name{getParmsLognormForMoments}
\alias{getParmsLognormForMoments}
\title{getParmsLognormForMoments}
\description{get the mean and variance of a log-normal distribution}
\usage{getParmsLognormForMoments(mean, var, sigmaOrig = sqrt(var))}
\arguments{
  \item{mean}{expected value at original scale}
  \item{var}{variance at original scale}
  \item{sigmaOrig}{standard deviation at original scale
, can be specified alternatively to the variance}
}

\value{numeric matrix with columns
\item{mu}{center parameter
(mean at log scale, log(median))}
\item{sigma}{scale parameter
(standard deviation at log scale)
}}
\references{\code{Limpert E, Stahel W & Abbt M (2001)
Log-normal Distributions across the Sciences: Keys and Clues.
Oxford University Press (OUP) 51, 341,
10.1641/0006-3568(2001)051[0341:lndats]2.0.co;2}}
\author{Thomas Wutzler}




\examples{
.mean <- 1
.var <- c(1.3,2)^2
getParmsLognormForMoments(.mean, .var)
}
