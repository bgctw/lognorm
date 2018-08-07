\name{getParmsLognormForMoments}
\alias{getParmsLognormForMoments}
\title{getParmsLognormForMoments}
\description{get the mean and variance of a log-normal distribution}
\usage{getParmsLognormForMoments(mean, var)}
\arguments{
  \item{mean}{expected value at original scale}
  \item{var}{variance at original scale}
}

\value{numeric matrix with columns
\item{mu}{center parameter
(mean at log scale, log(median))}
\item{sigma}{scale parameter
(standard deviation at log scale)
}}
\references{Limpert E, Stahel W & Abbt M (2001)
Log-normal Distributions across the Sciences: Keys and Clues.
Oxford University Press (OUP) 51, 341,
10.1641/0006-3568(2001)051[0341:lndats]2.0.co;2}
\author{Thomas Wutzler}





