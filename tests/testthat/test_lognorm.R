.tmp.f <- function(){
  require(testthat)  
}
context("lognormal")

test_that("getLognormMoments and getParmsLognormFromMoments vector",{
  mu = log(100)
  sigma = log(1.2)
  moments = getLognormMoments(mu, sigma)
  coefEst <- getParmsLognormForMoments( moments[,1], moments[,2])
  expect_equal( c(mu = mu, sigma = sigma), coefEst[1,])
})

test_that("getLognormMoments and getParmsLognormFromMoments matrix",{
  mu = log(c(100,200))
  sigma = rep(log(1.2),2)
  moments = getLognormMoments(mu, sigma)
  coefEst <- getParmsLognormForMoments( moments[,1], moments[,2])
  expect_equal( cbind(mu = mu, sigma = sigma), coefEst)
})

test_that("getParmsLognormFromExpValue matrix",{
  mu = log(c(100,200))
  sigma = rep(log(1.2),2)
  moments = getLognormMoments(mu, sigma)
  coefEst <- getParmsLognormForExpval( moments[,1], exp(sigma))
  expect_equal( cbind(mu = mu, sigma = sigma), coefEst)
})

test_that("getLognormMedian",{
  mu = log(c(100,200))
  sigma = rep(log(1.2),2)
  ans = getLognormMedian(mu, sigma)
  expect_equal( exp(mu), ans)
})

test_that("getLognormMode",{
  mu = log(c(100,200))
  sigma = rep(log(1.2),2)
  mode = getLognormMode(mu, sigma)
  expect_equal( exp(mu - sigma^2), mode)
})

test_that("scaleLogToOrig vector",{
  xLog <- data.frame(logmean = c(0.8, 0.8), sigma = c(0.2, 0.3))
  xOrig <- as.data.frame(scaleLogToOrig(xLog$logmean, xLog$sigma))
  xLog2 <- as.data.frame(scaleOrigToLog(xOrig$mean, xOrig$sd))
  expect_true( all.equal(xLog, xLog2) )
  xLog3 <- as.data.frame(getParmsLognormForMoments(xOrig$mean, xOrig$sd^2))
  expect_true( all.equal(xLog$sigma, xLog3$sigma) )
})

test_that("scaleLogToOrig scalar",{
  xLog <- data.frame(logmean = c(0.8), sigma = c(0.2))
  xOrig <- scaleLogToOrig(xLog$logmean, xLog$sigma)
  expect_equal(dim(xOrig), c(1,2))
  xLog2 <- scaleOrigToLog(xOrig[,"mean"], xOrig[,"sd"])
  expect_true( all.equal(xLog, as.data.frame(xLog2)) )
})




