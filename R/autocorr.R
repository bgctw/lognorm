
#' @export
computeEffectiveNumObs <- function(
  ### compute the effective number of observations taking into account autocorrelation
  res  ##<< numeric of autocorrelated numbers, usually observation - model residuals
  , effAcf = computeEffectiveAutoCorr(res) ##<< autocorrelation coefficients.
  ## The first entry is fixed at 1 for zero distance.
  ## May provide precomputed for efficiency or computed from a larger sample.
  , na.rm = FALSE  ##<< a logical value indicating whether NA values should be 
  ## stripped before the computation proceeds. 
){
  ##references<< 
  ## \code{Zieba & Ramza (2011) 
  ## Standard Deviation of the Mean of Autocorrelated 
  ## Observations Estimated with the Use of the Autocorrelation Function 
  ## Estimated From the Data. 
  ## Metrology and Measurement Systems, 
  ## Walter de Gruyter GmbH, 18 10.2478/v10178-011-0052-x}
  ## 
  ## \code{Bayley & Hammersley (1946) 
  ## The "effective" number of independent observations in an autocorrelated 
  ## time series. 
  ## Supplement to the Journal of the Royal Statistical Society, JSTOR,8,184-197}
  #
  ##details<< Handling of NA values: NAs at the beginning or end and are 
  ## just trimmed before computation and pose no problem. 
  ## However with NAs aside from edges, the return value is biased low,
  ## because correlation terms are subtracted for those positions.
  resTr <- .trimNA(res)
  if (!isTRUE(na.rm) & any(is.na(resTr))) return(NA_integer_)
  isFin <- is.finite(resTr)
  n <- sum(isFin)
  if (n == 0) return(0L)
  effAcfD <- effAcf[-1] # without zero distance
  # correlations may have been computed on a sample larger than given
  # then only use the 1..(n-1) components
  nC <- min(length(resTr) - 1, length(effAcfD))
  if (nC == 0) return(n)
  nEff0 <- n/(1 + 2*sum((1 - 1:nC/length(resTr))*effAcfD[1:nC]))
  ##details<< Because of NA correlation terms, the computed effective number of
  ## observations can be smaller than 1. In this case 1 is returned.
  nEff <- max(1, nEff0)  
  if ( nEff > n) stop("encountered nEff larger than finite records.")
  ##value<< integer scalar: effective number of observations
  nEff
}
attr(computeEffectiveNumObs,"ex") <- function(){
  # generate autocorrelated time series
  res <- stats::filter(rnorm(1000), filter = rep(1,5), circular = TRUE)
  res[100:120] <- NA
  # plot the series of autocorrelated random variables
  plot(res)
  # plot their empirical autocorrelation function
  acf(res, na.action = na.pass)
  #effAcf <- computeEffectiveAutoCorr(res)
  # the effective number of parameters is less than number of 1000 samples
  (nEff <- computeEffectiveNumObs(res, na.rm = TRUE))
}


.trimNA <- function(
  ### remove NA values at the start and end
  x  ##<< numeric vectpr
){
  nisna <- complete.cases(x)
  idx <- cumsum(nisna > 0) & rev(cumsum(rev(nisna))) > 0
  ##value<< subset of x with leading and trailing NAs removed
  x[idx]
}

#' @export
computeEffectiveAutoCorr <- function(
  ### return the vector of effective components of the autocorrelation
  res  ##<< numeric of autocorrelated numbers, usually observation - model residuals
){
  # first compute empirical autocorrelations
  ans <- acf(res, na.action = na.pass, plot = FALSE)
  # next get the number of elements before crossing the zero line
  nC <- suppressWarnings(min(which(ans$acf < 0)) - 1)
  if (!is.finite(nC)) {
    ans <- acf(res, na.action = na.pass, plot = FALSE, lag.max = Inf)
    # next get the number of elements before crossing the zero line
    nC <- suppressWarnings(min(which(ans$acf < 0)) - 1)
  }
  ##details<<
  ## Returns all components before first negative autocorrelation
  ##references<<
  ## \code{Zieba 2011 Standard Deviation of the Mean of Autocorrelated
  ## Observations Estimated with the Use of the Autocorrelation Function Estimated
  ## From the Data}
  n <- length(res)
  if (!is.finite(nC)) return(c(1))
  ##value<< numeric vector: strongest components of the autocorrelation function
  ans$acf[1:nC]
}
attr(computeEffectiveAutoCorr,"ex") <- function(){
  # generate autocorrelated time series
  res <- stats::filter(rnorm(1000), filter = rep(1,5), circular = TRUE)
  res[100:120] <- NA
  (effAcf <- computeEffectiveAutoCorr(res))
}


#' @export
varEffective <- function(
  ### estimate the variance of a correlated time series
  res  ##<< numeric of autocorrelated numbers, usually observation - model residuals
  , nEff = computeEffectiveNumObs(res, na.rm = na.rm) ##<< effective 
  ## number of observations
  , na.rm = FALSE  ##<< set to TRUE to remove NA cases before computation
  , ...  ##<< further arguments to \code{\link{var}}
) {
  ##details<< The BLUE is not anymore the usual variance, but a modified
  ## variance as given in \code{Zieba 2011}
  a <- 1
  var(res, na.rm = na.rm, ...) * nEff/(nEff - 1)
  ### The estimated variance of the sample
}
attr(varEffective,"ex") <- function(){
  # generate autocorrelated time series
  res <- stats::filter(rnorm(1000), filter = rep(1,5), circular = TRUE)
  res[100:120] <- NA
  # if correlations are neglected, the estimate of the variance is biased low
  (varNeglectCorr <- var(res, na.rm = TRUE))
  (varCorr <- varEffective(res, na.rm = TRUE))
}


.tmp.f <- function(){
  dss <- dsfP %>% mutate(
    isMeasured = (qfResp == 0)
    ,respPred = ifelse(isMeasured, computeRespLloydTaylor(
      dss$RRef, E0 = 308.56, temp = 273.15 + dss$temp, tempRef = 273.15 + 10), NA)
    , res = respPred - resp
  )
  plot( resp ~ date, filter(dss, isMeasured))
  points( respPred ~ date, filter(dss, isMeasured), col = "red")
  plot( res ~ date, dss)
  ans <- acf(dss$res, na.action = na.pass)
  nC <- min(which(ans$acf < 0)) - 1
  n <- nrow(dss)
  nEff <- n/(1 + 2*sum(1 - 1:nC/n*ans$acf[1:nC]))
  mean(dss$resp, na.rm = TRUE)
  sd(dss$resp, na.rm = TRUE)/sqrt(n)
  sd(dss$resp, na.rm = TRUE)/sqrt(nEff)
}

#' @export
getCorrMatFromAcf <- function(
  ### Construct the full correlation matrix from autocorrelation components.
  nRow      ##<< number of rows in correlation matrix
  , effAcf  ##<< numeric vector of effective autocorrelation components
  ##. The first entry, which is defined as 1, is not used.
){
  nDiag <- length(effAcf) - 1
  if (nDiag < 1) return(Diagonal(nRow))
  bandSparse(
    nRow, nRow, k = 0:nDiag, symmetric = TRUE
    , diagonals = lapply(0:nDiag, function(i){
      rep(effAcf[i + 1L], nRow - i)
    })
  )
}

#' @export
setMatrixOffDiagonals <- function(
  ### set off-diagonal values of the matrix
  x        ##<< numeric square matrix
  , diag = 1:length(value)  ##<< integer vector specifying the diagonals
  ## 0 is the center +1 the first
  ## row to upper and -2 the second row to lower
  , value ##<< numeric vector of values to fill in
  , isSymmetric = FALSE ##<< set to TRUE to to only
  ## specify the upper diagonal element but also
  ## set the lower in the mirrored diagonal
){
  if (!is.matrix(x) || !is.numeric(x) ) stop(
    "x must be a numeric matrix" )
  dimX <- dim(x)
  if (length(value) == 1) value <- rep(value, length(diag))
  if (length(value) != length(diag)) stop(
    "value must be of same length as diag")
  # A companion matrix that indicates how "off" a diagonal is:
  delta <- col(x) - row(x)
  if (isTRUE(isSymmetric)) {
    for (i in seq_along(diag)) {
      x[abs(delta) == diag[i] ] <- value[i]
    }
  } else {
    for (i in seq_along(diag)) {
      x[delta == diag[i] ] <- value[i]
    }
  }
  # dimensions may have gone lost
  dim(x) <- dimX
  ##value<< matrix with modified diagonal elements
  x
}
