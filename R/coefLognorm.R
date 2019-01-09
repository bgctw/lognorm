#' @export
getParmsLognormForMedianAndUpper <- function(
  ### Calculates mu and sigma of lognormal from median and upper quantile.
  median	##<< geometric mu (median at the original exponential scale)
  ,quant	##<< value at the upper quantile, i.e. practical maximum
  ,sigmaFac = qnorm(0.99) 	##<< sigmaFac=2 is 95% sigmaFac=2.6 is 99% interval
){
  # mu_geo = exp(mu); sigma_geo = exp(sigma)
  # logSpace cf-Interval: mu +- n sigma (n=1.96 for cf95)
  # expSpace cf-Interval: mu_geo */ simga_geo^n
  # given geometric mu (median at the original exponential scale)
  # and mu+sigmaFac*sigma
  tmp.mu = log(median)
  cbind(mu = tmp.mu, sigma = 1/sigmaFac*(log(quant) - tmp.mu) )
  ### named numeric vector: mu and sigma parameter of the lognormal distribution.
}

#' @export
getParmsLognormForMeanAndUpper <- function(
  ### Calculates mu and sigma of lognormal from median and upper quantile.
  mean	  ##<< expected value at the original scale
  ,quant	##<< value at the upper quantile, i.e. practical maximum
  ,sigmaFac = qnorm(0.99) 	##<< sigmaFac=2 is 95% sigmaFac=2.6 is 99% interval
){
  # solution of
  # (1) mle = exp(mu - sigma^2)
  # (2) quant = mu + sigmaFac sigma
  # see inst/doc/coefLognorm.Rmd for derivation
  ##details<< There are two valid solutions. This routine returns
  ## the one with lower sigma, i.e. the not so strongly skewed solution.
  m <- log(mean)
  #sigma <- sigmaFac + sqrt(sigmaFac^2 - 2*(log(quant) - m))
  sigma <- sigmaFac - sqrt(sigmaFac^2 - 2*(log(quant) - m))
  mu <- m - sigma^2/2
  cbind(mu = mu, sigma = sigma )
  ### numeric matrix: columns mu and sigma parameter of the
  ### lognormal distribution.
}



#' @export
getParmsLognormForLowerAndUpper <- function(
  ### Calculates mu and sigma of lognormal from lower and upper quantile.
  lower	  ##<< value at the lower quantile, i.e. practical minimum
  , upper	##<< value at the upper quantile, i.e. practical maximum
  , sigmaFac = qnorm(0.99) ##<< sigmaFac = 2 is 95%
  ## sigmaFac = 2.6 is 99% interval
  , isTransScale = FALSE ##<< if true lower and upper are already on log scale
){
  if (!isTRUE(isTransScale) ) {
    lower <- log(lower)
    upper <- log(upper)
  }
  sigma <- (upper - lower)/2/sigmaFac
  cbind( mu = upper - sigmaFac*sigma, sigma = sigma )
  ### named numeric vector: mu and sigma parameter of the lognormal
  ### distribution.
}
attr(getParmsLognormForLowerAndUpper,"ex") <- function(){
  # sample in normal space
  mu <- 5
  sigma <- 2
  rrNorm <- rnorm(1000, mean = mu, sd = sigma)
  # transform to orignal scale
  rrOrig <- exp(rrNorm)
  # and re-estimate parameters from original scale
  res <- getParmsLognormForMedianAndUpper(
    median(rrOrig), quantile(rrOrig, probs = 0.95), sigmaFac = qnorm(0.95))
  expected <- c(mu = mu,sigma = sigma)
  all.equal(res[1,], expected, tolerance = .1, scale = 1)
}


#' @export
getParmsLognormForModeAndUpper <- function(
  ### Calculates mu and sigma of lognormal from mode and upper quantile.
  mle			  ##<< numeric vector: mode at the original scale
  , quant		##<< numeric vector: value at the upper quantile,
  ## i.e. practical maximum
  , sigmaFac=qnorm(0.99) 	##<< sigmaFac=2 is 95% sigmaFac=2.6 is 99% interval
){
  # solution of
  # (1) mle = exp(mu - sigma^2)
  # (2) quant = mu + sigmaFac sigma
  # see inst/doc/coefLognorm.Rmd for derivation
  m <- log(mle)
  f2 <- sigmaFac/2
  sigma <- -f2 + sqrt(f2^2 + (log(quant) - m))
  mu <- m + sigma^2
  cbind(mu = mu, sigma = sigma )
  ### numeric matrix: columns mu and sigma parameter of the
  ### lognormal distribution.
  ### Rows correspond to rows of mle and quant
}
attr(getParmsLognormForModeAndUpper,"ex") <- function(){
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
  quant <- qlnorm(perc, meanlog = theta0[1], sdlog = theta0[2])
  (thetaEst <- getParmsLognormForModeAndUpper(mle,quant = quant,sigmaFac = qnorm(perc)) )

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

.ofLnormMLE <- function(
  ### Objective function used by \code{\link{coefLogitnormMLE}}.
  mu					  ##<< the mu parameter to estimate
  ,logMle				##<< the mode of the density distribution
  ,quant				##<< quantiles (values) at for percile perc
  ,perc = 0.975 ##<< percentiles where quant is given
){
  #mle = exp( mu-sigma^2)
  if (mu <= logMle) return(.Machine$double.xmax)
  sigma = sqrt(mu - logMle)
  tmp.predp = plnorm(quant, meanlog = mu, sdlog = sigma )
  tmp.diff = tmp.predp - perc
  sum(tmp.diff^2)
}
