#' Compute summary statistics of a log-normal distribution
#'
#' @describeIn getLognormMoments
#' get the expected value, variance, and coefficient of variation
#' @param mu numeric vector: location parameter
#' @param sigma numeric vector: scale parameter
#'
#' @return for \code{getLognormMoments} a numeric matrix with columns
#' \code{mean} (expected value at original scale)
#' , \code{var} (variance at original scale)
#' , and \code{cv} (coefficient of variation: sqrt(var)/mean).
#' For the other functions a numeric vector of the required summary.
#' 
#' @exampleFunction example_getLognormMoments
#' @references \code{Limpert E, Stahel W & Abbt M (2001)
#' Log-normal Distributions across the Sciences: Keys and Clues.
#' Oxford University Press (OUP) 51, 341,
#' 10.1641/0006-3568(2001)051[0341:lndats]2.0.co;2}
#' @export
getLognormMoments <- function(mu,sigma){
  sigma2 <- sigma*sigma
  m = exp(mu + sigma2/2)
  v = (exp(sigma2) - 1)*exp(2*mu + sigma2)
  cbind(
    mean = as.vector(m)              ##<< expected value at original scale
    , var = as.vector(v)             ##<< variance at original scale
    , cv = as.vector(sqrt(v)/m)      ##<< coefficient of variation: std/mean
  )
}
example_getLognormMoments <- function(){
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

#' @describeIn getLognormMoments
#' get the median
#' @export
getLognormMedian <- function(mu, sigma) exp(mu)

#' @describeIn getLognormMoments
#' get the mode
#' @export
getLognormMode <- function(mu,sigma) exp(mu - sigma*sigma)
