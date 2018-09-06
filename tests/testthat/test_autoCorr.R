.tmp.f <- function(){
  require(testthat)  
}
context("autoCorr")

test_that("getCorrMatFromAcf",{
  effAcf <- 1:5
  nRow <- 8
  m <- getCorrMatFromAcf(nRow, effAcf)
  expect_true(inherits(m,"Matrix"))
  expect_equal( dim(m), c(8,8))
  expect_equal( m[1,], c(effAcf, rep(0,3)))
  expect_equal( m[,2], c(effAcf[2], effAcf, rep(0,2)))
})

test_that("setMatrixOffDiagonals",{
  size <- 6
  mat <- diag(nrow = size, ncol = size)
  value = c(0.5,0.25,0.1)
  mat2 <- setMatrixOffDiagonals(
    mat, 1:3, value = value, isSymmetric =  TRUE)
  expect_true( all(c(mat2[1,2], mat2[2,1], mat2[3,2]) - value[1] == 0) )
  expect_true( all(c(mat2[1,3], mat2[3,1], mat2[4,2]) - value[2] == 0) )
  expect_true( all(c(mat2[1,5], mat2[5,1]) == 0) )
})

test_that(".trimNA",{
  x <- 1:5; x[2] <- NA
  expect_equal( .trimNA(x), x)
  y <- c(rep(NA,3),x )
  expect_equal( .trimNA(y), x)
  y <- c(x, rep(NA,3) )
  expect_equal( .trimNA(y), x)
  y <- c(rep(NA,3),x, rep(NA,3) )
  expect_equal( .trimNA(y), x)
})

test_that("computeEffectiveNumObs",{
  # generate autocorrelated time series
  res <- stats::filter(rnorm(1000), filter = rep(1,5), circular = TRUE)
  res[100:120] <- NA
  .tmp.f <- function(){
    plot(res)
    acf(res, na.action = na.pass)
  }
  effAcf <- computeEffectiveAutoCorr(res)
  #expect_true( length(effAcf) %in% 5:20) # depends on random numbers
  nEff <- computeEffectiveNumObs(res, na.rm = TRUE)
  expect_true(nEff < 1000 )
})

test_that("computeEffectiveNumObs with NA",{
  # generate autocorrelated time series
  res <- stats::filter(rnorm(1000), filter = rep(1,8), circular = TRUE)
  effAcfFull <- computeEffectiveAutoCorr(res)
  res[8:1000] <- NA
  nEff <- computeEffectiveNumObs(res)
  # if first correlation is negative, nEff will be equal
  expect_true(nEff <=  sum(is.finite(res))) 
  expect_true(nEff >=  1)
  res[3] <- NA
  nEff <- computeEffectiveNumObs(res)
  expect_true(is.na(nEff))
  nEff <- computeEffectiveNumObs(res, na.rm = TRUE)
  expect_true(nEff <=  sum(is.finite(res)))
  nEff <- computeEffectiveNumObs(res, na.rm = TRUE, effAcf = effAcfFull)
  expect_true(nEff <=  sum(is.finite(res)))
  expect_true(nEff >=  1)
  res[2:5] <- NA  # many NAs in center
  nEff <- computeEffectiveNumObs(res, na.rm = TRUE, effAcf = effAcfFull)
  expect_equal(nEff,1)
})

test_that("computeEffectiveNumObs with single finite obs",{
  # generate autocorrelated time series
  res <- stats::filter(rnorm(100), filter = rep(1,5), circular = TRUE)
  res[2:100] <- NA
  nEff <- computeEffectiveNumObs(res) # NA at edges is ok
  expect_equal(nEff,1)
  nEff <- computeEffectiveNumObs(res, na.rm = TRUE)
  expect_equal(nEff,1)
})

test_that("computeEffectiveNumObs with no finite obs",{
  # generate autocorrelated time series
  res <- stats::filter(rnorm(100), filter = rep(1,5), circular = TRUE)
  res[] <- NA
  nEff <- computeEffectiveNumObs(res)
  expect_equal(nEff,0)
  nEff <- computeEffectiveNumObs(res, na.rm = TRUE)
  expect_equal(nEff,0)
})

