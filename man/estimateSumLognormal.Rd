\name{estimateSumLognormal}
\alias{estimateSumLognormal}
\title{estimateSumLognormal}
\description{Estimate the distribution parameters of the lognormal approximation to the sum}
\usage{estimateSumLognormal(mu, sigma, corr = Diagonal(length(mu)), 
    sigmaSum = numeric(0), corrLength = if (inherits(corr, 
        "ddiMatrix")) 0 else nTerm, isStopOnNoTerm = FALSE, 
    effAcf, na.rm = isStopOnNoTerm)}
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
  \item{isStopOnNoTerm}{if no finite estimate is provided then by
default NA is returned for the sum.
Set this to TRUE to issue an error instead.}
  \item{effAcf}{numeric vector of effective autocorrelation
This overrides arguments \code{corr} and \code{corrLength}}
  \item{na.rm}{if there are terms with NA values in mu or sigma
by default also the sum coefficients are NA. Set to TRUE to 
neglect such terms in the sum.}
}

\value{numeric vector with two components mu and sigma
the parameters of the lognormal distribution at log scale}
\references{\code{Lo C (2013) WKB approximation for the sum of two correlated lognormal 
random variables.
Applied Mathematical Sciences, Hikari, Ltd., 7 , 6355-6367 
10.12988/ams.2013.39511}}
\author{Thomas Wutzler}




\examples{
# distribution of the sum of two lognormally distributed random variables
mu1 = log(110)
mu2 = log(100)
sigma1 = log(1.2)
sigma2 = log(1.6)
(coefSum <- estimateSumLognormal( c(mu1,mu2), c(sigma1,sigma2) ))
# repeat with correlation
(coefSumCor <- estimateSumLognormal( c(mu1,mu2), c(sigma1,sigma2), effAcf = c(1,0.9) ))
# expected value is equal, but variance with correlated variables is larger
getLognormMoments(coefSum["mu"],coefSum["sigma"])
getLognormMoments(coefSumCor["mu"],coefSumCor["sigma"])
}
