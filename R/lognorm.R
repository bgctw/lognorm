#' Compute summary statistics of a log-normal distribution
#'
#' @describeIn getLognormMoments
#' get the expected value, variance, and coefficient of variation
#'
#' @param mu numeric vector: location parameter
#' @param sigma numeric vector: scale parameter
#' @param m mean at original scale, may override default based on mu
#'
#' @return for \code{getLognormMoments} a numeric matrix with columns
#' \code{mean} (expected value at original scale)
#' , \code{var} (variance at original scale)
#' , and \code{cv} (coefficient of variation: sqrt(var)/mean).
#' For the other functions a numeric vector of the required summary.
#' 
#' @examples
#'   # start by estimating lognormal parameters from moments
#'   .mean <- 1
#'   .var <- c(1.3,2)^2
#'   parms <- getParmsLognormForMoments(.mean, .var)
#'   #
#'   # computed moments must equal previous ones
#'   (ans <- getLognormMoments(parms[,"mu"], parms[,"sigma"]))
#'   cbind(.var, ans[,"var"])
#'   #
#'   getLognormMedian(mu = log(1), sigma = log(2))
#'   getLognormMode(mu = log(1), sigma = c(log(1.2),log(2)))
#' 
#' @references \code{Limpert E, Stahel W & Abbt M (2001)
#' Log-normal Distributions across the Sciences: Keys and Clues.
#' Oxford University Press (OUP) 51, 341,
#' 10.1641/0006-3568(2001)051[0341:lndats]2.0.co;2}
#' @seealso scaleLogToOrig
#' @export
getLognormMoments <- function(mu,sigma, m = exp(mu + sigma2/2)){
  sigma2 <- sigma*sigma
  v = (exp(sigma2) - 1)*m^2
  cbind(
    mean = as.vector(m)              ##<< expected value at original scale
    , var = as.vector(v)             ##<< variance at original scale
    , cv = as.vector(sqrt(v)/m)      ##<< coefficient of variation: std/mean
  )
}

#' Scale standard deviation between log and original scale.
#'
#' When comparing values at log scale that have different sd at original scale, 
#' better compare log(mean) instead of mu.
#' 
#' @describeIn scaleLogToOrig get logmean and sigma at log scale
#' @param logmean log of the expected value
#' @param sigma standard deviation at log scale
#'
#' @return numeric matrix with columns \code{mean}, and \code{sd} at original scale
#' @examples
#'   xLog <- data.frame(logmean = c(0.8, 0.8), sigma = c(0.2, 0.3))
#'   xOrig <- as.data.frame(scaleLogToOrig(xLog$logmean, xLog$sigma))
#'   xLog2 <- as.data.frame(scaleOrigToLog(xOrig$mean, xOrig$sd))
#'   all.equal(xLog, xLog2)
#'   xLog3 <- as.data.frame(getParmsLognormForMoments(xOrig$mean, xOrig$sd^2))
#'   all.equal(xLog$sigma, xLog3$sigma) # but mu  < logmean 
#' @export
scaleLogToOrig <- function(logmean, sigma){
  ans <- moments <- getLognormMoments(sigma = sigma, 
                                      m = exp(logmean))[,1:2, drop = FALSE]
  colnames(ans) <- c("mean","sd")
  ans[,"sd"] <- sqrt(moments[,"var"])
  ans
}

#' @describeIn scaleLogToOrig get mean and sd at original scale
#' @param mean expected value at original scale
#' @param sd standard deviation at original scale
#' @export
scaleOrigToLog <- function(mean, sd){
  ans <- getParmsLognormForMoments(mean = mean, sigmaOrig = sd)
  colnames(ans) <- c("logmean","sigma")
  ans[,"logmean"] <- log(mean)
  ans
}

#' @describeIn getLognormMoments
#' get the median
#' @export
getLognormMedian <- function(mu, sigma) exp(mu)

#' @describeIn getLognormMoments
#' get the mode
#' @export
getLognormMode <- function(mu,sigma) exp(mu - sigma*sigma)
