#' @describeIn getParmsLognormForModeAndUpper 
#'   Calculates mu and sigma of lognormal from median and upper quantile.
#' @param median geometric mu (median at the original exponential scale)
#' @export
getParmsLognormForMedianAndUpper <- function(
  median, upper,sigmaFac = qnorm(0.99)){
  # mu_geo = exp(mu); sigma_geo = exp(sigma)
  # logSpace cf-Interval: mu +- n sigma (n=1.96 for cf95)
  # expSpace cf-Interval: mu_geo */ simga_geo^n
  # given geometric mu (median at the original exponential scale)
  # and mu+sigmaFac*sigma
  tmp.mu = log(median)
  cbind(mu = tmp.mu, sigma = 1/sigmaFac*(log(upper) - tmp.mu) )
}

#' @describeIn getParmsLognormForModeAndUpper
#'   Calculates mu and sigma of lognormal from mean and upper quantile.
#' @param mean expected value at the original scale
#' @details For \code{getParmsLognormForMeanAndUpper}
#' there are two valid solutions, and the one with lower sigma
#' , i.e. the not so strongly skewed solution is returned.
#' @export
getParmsLognormForMeanAndUpper <- function(
  mean, upper,sigmaFac = qnorm(0.99)
){
  # solution of
  # (1) mle = exp(mu - sigma^2)
  # (2) upper = mu + sigmaFac sigma
  # see inst/doc/coefLognorm.Rmd for derivation
  m <- log(mean)
  #sigma <- sigmaFac + sqrt(sigmaFac^2 - 2*(log(upper) - m))
  sigma <- sigmaFac - sqrt(sigmaFac^2 - 2*(log(upper) - m))
  mu <- m - sigma^2/2
  cbind(mu = mu, sigma = sigma )
}

#' @describeIn getParmsLognormForModeAndUpper
#'   Calculates mu and sigma of lognormal from lower and upper quantile.
#' @param lower value at the lower quantile, i.e. practical minimum
#' @param upper value at the upper quantile, i.e. practical maximum
#' @export
getParmsLognormForLowerAndUpper <- function(
  lower, upper, sigmaFac = qnorm(0.99)){
  getParmsLognormForLowerAndUpperLog(log(lower), log(upper), sigmaFac)
}

#' @describeIn getParmsLognormForModeAndUpper
#'   Calculates mu and sigma of lognormal from lower and upper quantile at log scale.
#' @param lowerLog value at the lower quantile, i.e. practical minimum at log scale
#' @param upperLog value at the upper quantile, i.e. practical maximum at log scale
#' @export
getParmsLognormForLowerAndUpperLog <- function(
  lowerLog, upperLog, sigmaFac = qnorm(0.99)){
  sigma <- (upperLog - lowerLog)/2/sigmaFac
  cbind(mu = upperLog - sigmaFac*sigma, sigma = sigma)
}

#' Calculate mu and sigma of lognormal from summary statistics.
#'
#' @describeIn getParmsLognormForModeAndUpper 
#'   Calculates mu and sigma of lognormal from mode and upper quantile.
#' @param mle numeric vector: mode at the original scale
#' @param upper numeric vector: value at the upper quantile,
#'    i.e. practical maximum
#' @param sigmaFac sigmaFac=2 is 95\% sigmaFac=2.6 is 99\% interval.
#'
#' @return numeric matrix wiht columns `mu` and `sigma`, the parameter of the
#'   lognormal distribution. Rows correspond to rows of inputs.
#' @export
#' @exampleFunction example_getParmsLognormForModeAndUpper
getParmsLognormForModeAndUpper <- function(
  mle, upper,sigmaFac = qnorm(0.99)){
  # solution of
  # (1) mle = exp(mu - sigma^2)
  # (2) upper = mu + sigmaFac sigma
  # see inst/doc/coefLognorm.Rmd for derivation
  m <- log(mle)
  f2 <- sigmaFac/2
  sigma <- -f2 + sqrt(f2^2 + (log(upper) - m))
  mu <- m + sigma^2
  cbind(mu = mu, sigma = sigma )
}
example_getParmsLognormForModeAndUpper <- function(){
  # example 1: a distribution with mode 1 and upper bound 5
  (thetaEst <- getParmsLognormForModeAndUpper(1,5))
  mle <- exp(thetaEst[1] - thetaEst[2]^2)
  all.equal(mle, 1, check.attributes = FALSE)

  # plot the distributions
  xGrid = seq(0,8, length.out = 81)[-1]
  dxEst <- dlnorm(xGrid, meanlog = thetaEst[1], sdlog = thetaEst[2])
  plot( dxEst~xGrid, type = "l",xlab = "x",ylab = "density")
  abline(v = c(1,5),col = "gray")

  # example 2: true parameters, which should be rediscovered
  theta0 <- c(mu = 1, sigma = 0.4)
  mle <- exp(theta0[1] - theta0[2]^2)
  perc <- 0.975		# some upper percentile, proxy for an upper bound
  upper <- qlnorm(perc, meanlog = theta0[1], sdlog = theta0[2])
  (thetaEst <- getParmsLognormForModeAndUpper(
    mle,upper = upper,sigmaFac = qnorm(perc)) )

  #plot the true and the rediscovered distributions
  xGrid = seq(0,10, length.out = 81)[-1]
  dx <- dlnorm(xGrid, meanlog = theta0[1], sdlog = theta0[2])
  dxEst <- dlnorm(xGrid, meanlog = thetaEst[1], sdlog = thetaEst[2])
  plot( dx~xGrid, type = "l")
  #plot( dx~xGrid, type = "n")
  #overplots the original, coincide
  lines( dxEst ~ xGrid, col = "red", lty = "dashed")

  # example 3: explore varying the uncertainty (the upper quantile)
  x <- seq(0.01,1.2,by = 0.01)
  mle = 0.2
  dx <- sapply(mle*2:8,function(q99){
    theta = getParmsLognormForModeAndUpper(mle,q99,qnorm(0.99))
    #dx <- dDistr(x,theta[,"mu"],theta[,"sigma"],trans = "lognorm")
    dx <- dlnorm(x,theta[,"mu"],theta[,"sigma"])
  })
  matplot(x,dx,type = "l")
}

#' @describeIn getParmsLognormForModeAndUpper 
#'   Calculate mu and sigma from moments (mean anc variance)
#'
#' @param mean expected value at original scale
#' @param var variance at original scale
#' @param sigmaOrig standard deviation at original scale 
#'
#' @references \code{Limpert E, Stahel W & Abbt M (2001)
#' Log-normal Distributions across the Sciences: Keys and Clues.
#' Oxford University Press (OUP) 51, 341,
#' 10.1641/0006-3568(2001)051[0341:lndats]2.0.co;2}
#' @export
getParmsLognormForMoments <- function(mean, var, sigmaOrig = sqrt(var)){
  omega = 1 + (sigmaOrig/mean)^2
  mu = log(mean / sqrt(omega))
  sigma = sqrt(log(omega))
  ans <- cbind(mu,sigma)
  rownames(ans) <- NULL # sometimes strange rownames by cbind
  ans
}


#' @describeIn getParmsLognormForModeAndUpper 
#'   Calculate mu and sigma from expected value and geometric standard deviation
#' @param median geometric mu (median at the original exponential scale)
#' @param sigmaStar multiplicative standard deviation
#' @exampleFunction example_getParmsLognormForExpval
#' @export
getParmsLognormForExpval <- function(mean, sigmaStar){
  sigma = log(sigmaStar)
  mu = log(mean) - sigma*sigma/2
  cbind(mu = mu, sigma = sigma)
}
example_getParmsLognormForExpval <- function(){
  # Calculate mu and sigma from expected value and geometric standard deviation
  .mean <- 1
  .sigmaStar <- c(1.3,2)
  (parms <- getParmsLognormForExpval(.mean, .sigmaStar))
  # multiplicative standard deviation must equal the specified value
  cbind(exp(parms[,"sigma"]), .sigmaStar)
}


#' Estimate lognormal distribution parameters from a sample
#'
#' @describeIn estimateParmsLognormFromSample 
#'    Estimate lognormal distribution parameters from a sample
#'
#' @param x numeric vector of sampled values
#' @param na.rm a logical value indicating whether 
#' NA values should be stripped before the computation proceeds.
#'
#' @return numeric vector with components \code{mu} and \code{sigma},
#' ie.., the center parameter (mean at log scale, log(median)) and 
#' the scale parameter (standard deviation at log scale)
#' 
#' @examples 
#' .mu <- log(1)
#' .sigma <- log(2)
#' n = 200
#' x <- exp(rnorm(n, mean = .mu, sd = .sigma))
#' exp(pL <- estimateParmsLognormFromSample(x)) # median and multiplicative stddev
#' c(mean(x), meanx <- getLognormMoments(pL["mu"],pL["sigma"])[,"mean"])
#' c(sd(x), sdx <- sqrt(getLognormMoments(pL["mu"],pL["sigma"])[,"var"]))
#' 
#' # stddev decreases (each sample about 0.9) to about 0.07
#' # for the mean with n replicated samples
#' se <- estimateStdErrParms(x)
#' sqrt(getLognormMoments(se["mu"],se["sigma"])[,"var"])
#' sd(x)/sqrt(n-1) # well approximated by normal
#' # expected value stays the same
#' c(meanx, getLognormMoments(se["mu"],se["sigma"])[,"mean"])
#' @export
estimateParmsLognormFromSample <- function(x, na.rm = FALSE){
  logx <- log(x)
  c(mu = mean(logx, na.rm = na.rm) ,sigma = sd(logx, na.rm = na.rm))
}

#' @describeIn estimateParmsLognormFromSample 
#'    Estimate parameters of the lognormal distribution of the mean from an uncorrelated sample
#' @details The expected value of a can be determined with
#'   higher accuracy the larger the sample. Here, the uncorrelated
#'   assumption is applied at the log scale and distribution parameters
#'   are returned with the same expected value as the sample, but with
#'   uncertainty (sigma) decreased by sqrt(nfin - 1). 
#'   
#'   Since with low relative error, the lognormal becomes very close
#'   to the normal distribution, the distribution of the mean can be
#'   well approximated by a normal with sd(mean(x)) ~ sd(x)/sqrt(n-1).
#' 
#' @export
estimateStdErrParms <- function(x, na.rm = FALSE){
  pl <- estimateParmsLognormFromSample(x = x, na.rm = na.rm)
  sigma <- pl["sigma"]
  logmean <- pl["mu"] + sigma^2/2
  nfin <- sum(is.finite(x))
  sdmean <- sigma / sqrt(nfin - 1)
  setNames(c(logmean - sdmean^2/2, sdmean),c("mu","sigma"))
}



