.tmp.f <- function(){
  require(testthat)
  #
  require(Matrix)
  require(tidyr)
  require(dplyr)
  require(ggplot2)
}
context("diffLognormal")

test_that("estimateDiffLognormal two Vars",{
  # generate nSample values of two lognormal random variables
  mu1 = log(110)
  mu2 = log(100)
  sigma1 = 0.25
  sigma2 = 0.15
  #estimateSumLognormalBenchmark( c(mu1,mu2), c(sigma1,sigma2) )
  coefSum <- estimateDiffLognormal( mu1,mu2,sigma1,sigma2 )
  m <- getLognormMoments(coefSum["mu"], coefSum["sigma"], 
                        shift = coefSum["shift"])[,"mean"]
  expect_equal(m, c(mean = exp(mu1 + sigma1^2/2) - exp(mu2 + sigma2^2/2)))
  # regression test
  coefSumExp <- c(mu = 6.15, sigma = 0.069, shift = 456.07)
  expect_equal( coefSum, coefSumExp, tolerance = 0.02 )
})

.tmp.plotsample <- function(){
  nSample = 500
  #nSample = 1e5
  ds <- data.frame(
    x1 = rlnorm(nSample, mu1, sigma1)
    , x2 = rlnorm(nSample, mu2, sigma2)
  ) %>%  mutate(
    y = x1 - x2
  )
  c(mean(ds$y), sd(ds$y))
  dsw <- gather(ds, "var", "value", x1, x2, y)
  p1 <- ggplot(dsw, aes(value, color = var)) + geom_density(linetype = "dotted")
  #
  p <- seq(0,1,length.out = 32)[-c(1,32)]
  dsPred1 <- data.frame(
    var = "x1", q = qlnorm(p, mu1, sigma1 )
  ) %>%
    mutate( d = dlnorm(q, mu1, sigma1))
  dsPred2 <- data.frame(
    var = "x2", q = qlnorm(p, mu2, sigma2 )
  ) %>%
    mutate( d = dlnorm(q, mu2, sigma2))
  dsPred <- rbind(dsPred1, dsPred2)
  p1 + geom_line(data = dsPred, aes(q, d))
  #
  coefSum <- estimateDiffLognormal( mu1,mu2,sigma1,sigma2 )
  dsPredY <- data.frame(
    var = "y", 
    q_shifted = qlnorm(p, coefSum["mu"], coefSum["sigma"] )
  ) %>%
    mutate( 
      q = q_shifted - coefSum["shift"],
      d = dlnorm(q_shifted, coefSum["mu"], coefSum["sigma"])
    )
  dsPred <- rbind(dsPred1, dsPred2, select(dsPredY, var, q, d))
  p1 + geom_line(data = dsPred, aes(q, d))
}


test_that("estimateDiffLognormal two Vars positive correlation",{
  # generate nSample values of two lognormal random variables
  mu1 = log(110)
  mu2 = log(100)
  sigma1 = 0.25
  sigma2 = 0.15
  #estimateSumLognormalBenchmark( c(mu1,mu2), c(sigma1,sigma2) )
  coefSum <- estimateDiffLognormal( mu1,mu2,sigma1,sigma2, corr = 0.8 )
  coefSumExp <- c(mu = 4.99, sigma = 0.13, shift = 134.14 )
  expect_equal( coefSum, coefSumExp, tolerance = 0.02 )
})

test_that("estimateDiffLognormal two Vars negative correlation",{
  # generate nSample values of two lognormal random variables
  mu1 = log(110)
  mu2 = log(100)
  sigma1 = 0.25
  sigma2 = 0.15
  #estimateSumLognormalBenchmark( c(mu1,mu2), c(sigma1,sigma2) )
  coefSum <- estimateDiffLognormal( mu1,mu2,sigma1,sigma2, corr = -0.8 )
  coefSumExp <- c(mu = 6.67, sigma = 0.053, shift = 778.01)
  expect_equal( coefSum, coefSumExp, tolerance = 0.02 )
})

test_that("pDiffLognormalSample",{
  # generate nSample values of two lognormal random variables
  mu1 = log(120)
  mu2 = log(60)
  sigma1 = 0.25
  sigma2 = 0.15
  #estimateSumLognormalBenchmark( c(mu1,mu2), c(sigma1,sigma2) )
  coefSum <- estimateDiffLognormal( mu1,mu2,sigma1,sigma2, corr = -0.8 )
  coefSum
  q <- c(-1000,0,2,1000)
  pLo <- plnorm(q + coefSum["shift"], coefSum["mu"], coefSum["sigma"])
  #.tmp.plotsample()
  pSample <- pDiffLognormalSample(mu1,mu2,sigma1,sigma2, corr = -0.8, q = q)
  #cbind(pLo, pSample)
  expect_equal( pLo, pSample, tolerance = 0.05)
})

