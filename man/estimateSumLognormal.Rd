\name{estimateSumLognormal}
\alias{estimateSumLognormal}
\title{estimateSumLognormal}
\description{Estimate the distribution parameters of the lognormal approximation to the sum}
\usage{estimateSumLognormal(mu, sigma, corr = Diagonal(length(mu)), 
    sigmaSum = numeric(0), corrLength = if (inherits(corr, 
        "ddiMatrix")) 0 else nTerm, isStopOnNoTerm = FALSE, 
    effAcf)}
\arguments{
  \item{mu}{numeric vector of center parameters of terms at log scale}
  \item{sigma}{numeric vector of variance parameter of terms at log scale}
  \item{corr}{numeric matrix
of correlations between the random variables}
  \item{sigmaSum}{numeric scalar: possibility to specify
of a precomputed scale parameter}
  \item{corrLength}{integer
scalar: set correlation length to smaller values
to speed up computation by neglecting correlations among terms
further apart.
Set to zero to omit correlations.}
  \item{isStopOnNoTerm}{if no finite estiamte is provide, by
default, NA is returned for the sum.
Set this to TRUE to issue an error instead.}
  \item{effAcf}{numeric vector of effective autocorrelation
This overides arguments \code{corr} and \code{corrLength}}
}

\value{numeric vector with two components mu and sigma
the parameters of the lognormal distribution at log scale}
\references{Lo C (2013) WKB approximation for the sum of two correlated lognormal 
random variables.
Applied Mathematical Sciences, Hikari, Ltd., 7 , 6355-6367 
10.12988/ams.2013.39511}
\author{Thomas Wutzler}





