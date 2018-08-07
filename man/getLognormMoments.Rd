\name{getLognormMoments}
\alias{getLognormMoments}
\title{getLognormMoments}
\description{get the expected value and variance of a log-normal distribution}
\usage{getLognormMoments(mu, sigma)}
\arguments{
  \item{mu}{center parameter (mean at log scale, log(median))}
  \item{sigma}{scale parameter (standard deviation at log scale)}
}

\value{numeric matrix with columns
\item{mean}{expected value at original scale}
\item{var}{variance at original scale}
\item{cv}{coefficient of variation: std/mean
}}
\references{Limpert E, Stahel W & Abbt M (2001)
Log-normal Distributions across the Sciences: Keys and Clues.
Oxford University Press (OUP) 51, 341,
10.1641/0006-3568(2001)051[0341:lndats]2.0.co;2}
\author{Thomas Wutzler}





