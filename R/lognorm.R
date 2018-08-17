#' @export
getLognormMoments <- function(
  ### get the expected value and variance of a log-normal distribution
  mu        ##<< center parameter (mean at log scale, log(median))
  , sigma   ##<< scale parameter (standard deviation at log scale)
){
  ##references<< Limpert E, Stahel W & Abbt M (2001)
  ## Log-normal Distributions across the Sciences: Keys and Clues.
  ## Oxford University Press (OUP) 51, 341,
  ## 10.1641/0006-3568(2001)051[0341:lndats]2.0.co;2
  sigma2 <- sigma*sigma
  m = exp(mu + sigma2/2)
  v = (exp(sigma2) - 1)*exp(2*mu + sigma2)
  ##value<< numeric matrix with columns
  cbind(
    mean = as.vector(m)              ##<< expected value at original scale
    , var = as.vector(v)             ##<< variance at original scale
    , cv = as.vector(sqrt(v)/m)      ##<< coefficient of variation: std/mean
  )
}

#' @export
getLognormMedian <- function(
  ### get the median of a log-normal distribution
  mu        ##<< center parameter (mean at log scale, log(median))
  , sigma = NA ##<< dummy not used, but signature as with Mode and moments
){
  ##value<< the median
  exp(mu)
}

#' @export
getLognormMode <- function(
  ### get the mode of a log-normal distribution
  mu        ##<< center parameter (mean at log scale, log(median))
  , sigma   ##<< scale parameter (standard deviation at log scale)
){
  ##value<< the mode
  exp(mu - sigma*sigma)
}

#' @export
getParmsLognormForMoments <- function(
  ### get the mean and variance of a log-normal distribution
  mean        ##<< expected value at original scale
  , var       ##<< variance at original scale
  , sigmaOrig = sqrt(var)  ##<< alternatively to the variance,
  ## the standard devation at original scale can be given
){
  ##references<< Limpert E, Stahel W & Abbt M (2001)
  ## Log-normal Distributions across the Sciences: Keys and Clues.
  ## Oxford University Press (OUP) 51, 341,
  ## 10.1641/0006-3568(2001)051[0341:lndats]2.0.co;2
  omega = 1 + (sigmaOrig/mean)^2
  mu = log(mean / sqrt(omega))
  sigma = sqrt(log(omega))
  ##value<< numeric matrix with columns
  cbind(
    mu = mu         ##<< center parameter
    ## (mean at log scale, log(median))
    ,sigma = sigma  ##<< scale parameter
    ## (standard deviation at log scale)
  )
}

#' @export
getParmsLognormForExpval <- function(
  ### get the lognormal parameters by expected value
  mean        ##<< expected value at original scale
  , sigmaStar ##<< multiplicative standard deviation
){
  sigma = log(sigmaStar)
  mu = log(mean) - sigma*sigma/2
  cbind(
    mu = mu         ##<< center parameter
    ## (mean at log scale, log(median))
    ,sigma = sigma  ##<< scale parameter
    ## (standard deviation at log scale)
  )
}

#' @export
estimateParmsLognormFromSample <- function(
  ### get the lognormal parameters by expected value.
  x        ##<< numeric vector of sampled values
  , na.rm = FALSE ##<< a logical value indicating whether 
  ## NA values should be stripped before the computation proceeds.
){
  logx <- log(x)
  c(
    mu = mean(logx, na.rm = na.rm)         ##<< center parameter
    ## (mean at log scale, log(median))
    ,sigma = sd(logx, na.rm = na.rm)  ##<< scale parameter
    ## (standard deviation at log scale)
  )
}
