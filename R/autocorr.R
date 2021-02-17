#' Compute the standard error accounting for empirical autocorrelations
#'
#' @param x numeric vector
#' @param na.rm logical. Should missing values be removed?
#' @param effCov numeric vector of effective covariance components
#'  first entry is the variance. See \code{\link{computeEffectiveAutoCorr}}
#' @details Computation follows 
#'  https://stats.stackexchange.com/questions/274635/calculating-error-of-mean-of-time-series.
#'  
#' @details The default uses empirical autocorrelation
#'  estimates from the supplied data up to first negative component.
#'  For short series of \code{x} it is strongly recommended to to
#'  provide \code{effCov} that was estimated on a longer time series.
#'
#' @export
#' @return numeric scalar of standard error of the mean of x
seCor <- function(
  x  
  , na.rm = FALSE 
  , effCov = computeEffectiveAutoCorr(x, type = "covariance")
){
  n <- if (na.rm) length(na.omit(x)) else length(x)
  #effAcf <- computeEffectiveAutoCorr(x, na.action = {if(na.rm) na.omit else na.pass})
  # do not remove NAs for autocorrelation computation to preserve distances
  if (n == 0) return(NA_real_)
  g1 <- effCov[1:min(length(x),length(effCov))]
  kmax <- length(g1) - 1
  # if there is no empirical autocorrelation
  if (kmax == 0) {
    varx <- var(x, na.rm = na.rm)
    return(sqrt(varx/n))
  } 
  g0 <- g1[1]
  g <- g1[-1]
  k <- 1:kmax
  varCor <- 1/n*(g0 + 2*sum( (n - k)/n * g))
  sqrt(varCor)
}

#' Compute the effective number of observations taking into account autocorrelation
#' 
#' @param res numeric of autocorrelated numbers, usually observation -
#'  model residuals
#' @param effAcf autocorrelation coefficients.
#'  The first entry is fixed at 1 for zero distance.
#' @param na.rm a logical value indicating whether NA values should be 
#'  stripped before the computation proceeds.
#'
#' @references 
#' \code{Zieba & Ramza (2011) 
#' Standard Deviation of the Mean of Autocorrelated 
#' Observations Estimated with the Use of the Autocorrelation Function 
#' Estimated From the Data. 
#' Metrology and Measurement Systems, 
#' Walter de Gruyter GmbH, 18 10.2478/v10178-011-0052-x}
#'
#' \code{Bayley & Hammersley (1946) 
#' The "effective" number of independent observations in an autocorrelated 
#' time series. 
#' Supplement to the Journal of the Royal Statistical Society, JSTOR,8,184-197}
#'
#' @details Handling of NA values: NAs at the beginning or end are 
#' just trimmed before computation and pose no problem. 
#' However with NAs aside from edges, the return value of nEff is biased low,
#' because correlation for pairs involving NA is still accounted in the 
#' denominator of (3) in Zieba 2011.
#' This leads to a conservative (biased high) estimates of standard errors.
#' The effect should be small if length(acf) < n_finite.
#' @details Because of NA correlation terms, the computed effective number of
#' observations can be smaller than 1. In this case 1 is returned.

#' @export
#' @return integer scalar: effective number of observations
#' @exampleFunction example_computeEffectiveNumObs
computeEffectiveNumObs <- function(
  res  
  , effAcf = computeEffectiveAutoCorr(res)
  , na.rm = FALSE  
){
  resTr <- .trimNA(res)
  if (!isTRUE(na.rm) & any(is.na(resTr))) return(NA_integer_)
  isFin <- is.finite(resTr)
  n <- sum(isFin)
  if (n == 0) return(0L)
  effAcfD <- effAcf[-1] # without zero distance
  # correlations may have been computed on a sample larger than given
  # then only use the 1..(n-1) components
  # If there an NAs in resTr then nFin < length(resTr) and the numerator is 
  # decreased.
  # But autocorrelation is still accounted for pairs involving NA.
  # Therefore, the denominator is overestimated and hence the fraction nEff is 
  # underestimated.
  nC <- min(length(resTr) - 1, length(effAcfD))
  #?nC <- min(n - 1, length(effAcfD))
  if (nC == 0) return(n)
  nEff0 <- n/(1 + 2*sum((1 - 1:nC/length(resTr))*effAcfD[1:nC]))
  # alternative for missings: using n=nFin instead of length(resTr) in k/n
  # i.e. computed on vector omitting NA
  # k/nFin is smaller then k/length(resTr) -> (1-k/N) larger -> 
  #   denominator larger -> nEff smaller 
  # this increases the bias underestimating nEff -> prefer length(nTr) solution
  nEff <- max(1, nEff0)  
  if ( nEff > n) stop("encountered nEff larger than finite records.")
  nEff
}
example_computeEffectiveNumObs <- function(){
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

#' Estimate vector of effective components of the autocorrelation
#' 
#' @param res  numeric of autocorrelated numbers, usually observation - 
#'  model residuals
#' @param type type of residuals (see \code{\link{acf}})
#'
#' @details
#'  Returns all components before first negative autocorrelation
#' @references
#'  \code{Zieba 2011 Standard Deviation of the Mean of Autocorrelated
#'  Observations Estimated with the Use of the Autocorrelation Function Estimated
#'  From the Data}
#' @export
#' @return numeric vector: strongest components of the autocorrelation function
#' @exampleFunction example_computeEffectiveAutoCorr
computeEffectiveAutoCorr <- function(
  res, type = "correlation" ){
  # first compute empirical autocorrelations
  ans <- acf(res, na.action = na.pass, plot = FALSE, type = type)
  # next get the number of elements before crossing the zero line
  nC <- suppressWarnings(min(which(ans$acf <= 0)) - 1)
  if (!is.finite(nC)) {
    # if there was no below zero correlation within defalt lag.max then
    # repeat acf with compting all lags
    ans <- acf(res, na.action = na.pass, plot = FALSE, type = type, lag.max = Inf)
    # append -1 so that nC equals to full lenght if no negative correlation
    nC <- min(which(c(ans$acf,-1) <= 0)) - 1
  }
  nC <- pmax(1,nC) # return at least one component
  ans$acf[1:nC]
}
example_computeEffectiveAutoCorr <- function(){
  # generate autocorrelated time series
  res <- stats::filter(rnorm(1000), filter = rep(1,5), circular = TRUE)
  res[100:120] <- NA
  (effAcf <- computeEffectiveAutoCorr(res))
}




#' remove NA values at the start and end
#'
#' @param x numeric vector
#'
#' @return subset of x with leading and trailing NAs removed
#' @keywords internal
.trimNA <- function(x){
  nisna <- complete.cases(x)
  idx <- cumsum(nisna > 0) & rev(cumsum(rev(nisna))) > 0
  x[idx]
}


#' Estimate the variance of a correlated time series
#' 
#' @param res numeric of autocorrelated numbers, usually observation - 
#'  model residuals
#' @param nEff effective number of observations, can be specified for efficiency
#' @param na.rm set to TRUE to remove NA cases before computation
#' @param ... urther arguments to \code{\link{var}}
#'
#' @export
#' @exampleFunction example_varEffective
#' @details The BLUE is not anymore the usual variance, but a modified
#'  variance as given in \code{Zieba 2011}
#'  @return The estimated variance of the sample
varEffective <- function(
  res
  , nEff = computeEffectiveNumObs(res, na.rm = na.rm) 
  , na.rm = FALSE  
  , ...  
) {
  var(res, na.rm = na.rm, ...) * nEff/(nEff - 1)
}
example_varEffective <- function(){
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

#' Construct the full correlation matrix from autocorrelation components.
#'
#' @param nRow number of rows in correlation matrix
#' @param effAcf numeric vector of effective autocorrelation components
#' . The first entry, which is defined as 1, is not used.
#'
#' @export
getCorrMatFromAcf <- function(nRow, effAcf){
  nDiag <- length(effAcf) - 1
  if (nDiag < 1) return(Diagonal(nRow))
  bandSparse(
    nRow, nRow, k = 0:nDiag, symmetric = TRUE
    , diagonals = lapply(0:nDiag, function(i){
      rep(effAcf[i + 1L], nRow - i)
    })
  )
}

#' set off-diagonal values of a matrix
#'
#' @param x numeric square matrix
#' @param diag integer vector specifying the diagonals
#'  0 is the center +1 the first row to upper and -2 the second row to lower
#' @param value numeric vector of values to fill in
#' @param isSymmetric  set to TRUE to to only
#'  specify the upper diagonal element but also
#'  set the lower in the mirrored diagonal
#'
#' @export
#' @return matrix with modified diagonal elements
setMatrixOffDiagonals <- function(
  x, diag = 1:length(value), value, isSymmetric = FALSE){
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
  x
}

