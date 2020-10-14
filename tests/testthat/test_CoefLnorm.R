.tmp.f <- function(){
  require(testthat)
}

context("getParmsLognormForMedianAndUpper")

test_that("getParmsLognormForMedianAndUpper",{
  mu = 2
  sd = c(1,0.8)
  p = 0.99
  med <- qlnorm(0.5, mu, sd)
  upper <- u <- qlnorm(p, mu, sd )		# p-confidence interval
  cf <- getParmsLognormForMedianAndUpper(med, upper)
  expect_equal( cf[,"mu"] , c(mu,mu) )
  expect_equal( cf[,"sigma"] , sd )
})

test_that("getParmsLognormForLowerAndUpper",{
  mu = 2
  sd = c(1,0.8)
  p = 0.99
  lower <- l <- qlnorm(1 - p, mu, sd )		# p-confidence interval
  upper <- u <- qlnorm(p, mu, sd )		# p-confidence interval
  cf <- getParmsLognormForLowerAndUpper(lower,upper)
  expect_equal( cf[,"mu"] , c(mu,mu) )
  expect_equal( cf[,"sigma"] , sd )
})

test_that("getParmsLognormForModeAndUpper",{
  thetaEst <- getParmsLognormForModeAndUpper(1,5)
  mle <- exp(thetaEst[1] - thetaEst[2]^2)
  expect_equal(mle , 1, check.attributes = FALSE)
  #
  q <- c(2,7)
  #trace(coefLognormMLE, recover)     #untrace(coefLognormMLE)
  theta <- getParmsLognormForModeAndUpper(q[1],q[2])
  tmp <- plnorm(q[2], meanlog = theta[1], sdlog = theta[2])
  q2 <- qlnorm(c(0.5,0.99), meanlog = theta[1], sdlog = theta[2])
  expect_equal(q[2],q2[2])
  mle <- exp(theta[1] - theta[2]^2)
  expect_equal(q[1],mle)
  #
  res <- getParmsLognormForModeAndUpper(1:5,q[2])
  expect_equal(5,nrow(res))
  expect_equal(res[2,,drop = FALSE], theta)
  res <- getParmsLognormForModeAndUpper(q[1],q[2] + (-2:2))
  expect_equal(5,nrow(res))
  expect_equal(res[3,,drop = FALSE], theta)
})

test_that("getParmsLognormForMeanAndUpper",{
  #trace(getParmsLognormForMeanAndUpper, recover) # untrace(getParmsLognormForMeanAndUpper)
  thetaEst <- getParmsLognormForMeanAndUpper(1,5)
  mean <- exp(thetaEst[1] + thetaEst[2]^2/2)
  expect_equal(mean , 1, check.attributes = FALSE)
  #
  q <- c(2,7)
  #trace(coefLognormMLE, recover)     #untrace(coefLognormMLE)
  theta <- getParmsLognormForMeanAndUpper(q[1],q[2])
  q2 <- qlnorm(c(0.5,0.99), meanlog = theta[1], sdlog = theta[2])
  expect_equal(q[2],q2[2])
  mean <- exp(theta[1] + theta[2]^2/2)
  expect_equal(q[1],mean)
  #
  res <- getParmsLognormForMeanAndUpper(1:5,q[2])
  expect_equal(5,nrow(res))
  expect_equal(res[2,,drop = FALSE], theta)
  res <- getParmsLognormForMeanAndUpper(q[1],q[2] + (-2:2))
  expect_equal(5,nrow(res))
  expect_equal(res[3,,drop = FALSE], theta)
})

test_that("estimateStdErrParms",{
  meanx <- 2
  relerr <- c(0.5)
  pS <- getParmsLognormForMoments(meanx, sigmaOrig = meanx*relerr)
  exp(pS)
  n <- 512
  n <- 32
  nrep <- 4048*64
  nrep <- 4048
  # bootstrap sample of mean across n lognormal variables
  mx_boot = purrr::map_dbl(seq_len(nrep), function(i){
    x <- exp(rnorm(n, mean = pS[1,"mu"], sd = pS[1,"sigma"]))
    #plot(density(x))
    #estimateParmsLognormFromSample(x)
    mean(x)
  })
  sdm <- sd(mx_boot)
  x <- exp(rnorm(n, mean = pS[1,"mu"], sd = pS[1,"sigma"]))
  ans <- estimateStdErrParms(x)
  xo <- scaleLogToOrig(ans["mu"], ans["sigma"])
  .tmp.f <- function(){
    exp(ans)
    xbounds <- exp(qnorm(c(0.025, 0.975), ans["mu"], ans["sigma"]))
    xgrid <- seq(xbounds[1]*0.9, xbounds[2]*1.1, length.out = 300)
    # orange: density of a deviation from true mean, corrected by 
    # difference between original mean and sample mean (meanx - mean(x))
    dans <- dnorm(log(xgrid - (meanx - mean(x))), ans["mu"], ans["sigma"])/mean(x)
    # blue: assuming normal distribution of x (also corrected by the sample mean)
    dansN <- dnorm(xgrid - (meanx - mean(x)), mean(x), sd(x)/sqrt(n - 1))
    plot(density(mx_boot)); 
    lines(dans ~ xgrid, col = "blue")
    lines(dansN ~ xgrid, col = "orange"); abline(v = meanx)
    # orange slightly better approximation, but here normal approx is sufficient
    # (with n and nrep sufficiently large)
  }
  expect_equivalent(mean(x), xo[1,"mean"], 0.2/sqrt(n), scale = 1 )
  expect_equivalent(sdm, xo[1,"sd"], 0.5/sqrt(n), scale = 1)
  # 
})








