% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lognorm.R
\name{getLognormMoments}
\alias{getLognormMoments}
\alias{getLognormMedian}
\alias{getLognormMode}
\title{Compute summary statistics of a log-normal distribution}
\usage{
getLognormMoments(mu, sigma, m = exp(mu + sigma2/2) - shift, shift = 0)

getLognormMedian(mu, sigma, shift = 0)

getLognormMode(mu, sigma, shift = 0)
}
\arguments{
\item{mu}{numeric vector: location parameter}

\item{sigma}{numeric vector: scale parameter}

\item{m}{mean at original scale, may override default based on mu}

\item{shift}{shift for the shifted lognormal distribution}
}
\value{
for \code{getLognormMoments} a numeric matrix with columns
\code{mean} (expected value at original scale)
, \code{var} (variance at original scale)
, and \code{cv} (coefficient of variation: sqrt(var)/mean).
For the other functions a numeric vector of the required summary.
}
\description{
Compute summary statistics of a log-normal distribution
}
\section{Functions}{
\itemize{
\item \code{getLognormMoments}: get the expected value, variance, and coefficient of variation

\item \code{getLognormMedian}: get the median

\item \code{getLognormMode}: get the mode
}}

\examples{
  # start by estimating lognormal parameters from moments
  .mean <- 1
  .var <- c(1.3,2)^2
  parms <- getParmsLognormForMoments(.mean, .var)
  #
  # computed moments must equal previous ones
  (ans <- getLognormMoments(parms[,"mu"], parms[,"sigma"]))
  cbind(.var, ans[,"var"])
  #
  getLognormMedian(mu = log(1), sigma = log(2))
  getLognormMode(mu = log(1), sigma = c(log(1.2),log(2)))

}
\references{
\code{Limpert E, Stahel W & Abbt M (2001)
Log-normal Distributions across the Sciences: Keys and Clues.
Oxford University Press (OUP) 51, 341,
10.1641/0006-3568(2001)051[0341:lndats]2.0.co;2}
}
\seealso{
scaleLogToOrig
}
