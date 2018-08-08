.tmp.f <- function(){
  require(testthat)  
}
context("autoCorr")

test_that("computeEffectiveNumObs",{
  # generate autocorrelated time series
  res <- stats::filter(rnorm(1000), filter = rep(1,5), circular = TRUE)
  res[100:120] <- NA
  .tmp.f <- function(){
    plot(res)
    acf(res, na.action = na.pass)
  }
  effAcf <- computeEffectiveAutoCorr(res)
  expect_true( length(effAcf) %in% 5:20) # depends on random numbers
  nEff <- computeEffectiveNumObs(res)
  expect_true(nEff < 1000 )
})

