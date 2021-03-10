#' Compute the standard error accounting for empirical autocorrelations
#'
#' @param x numeric vector
#' @param na.rm logical. Should missing values be removed?
#' @param effCor numeric vector of effective correlation components
#'  first entry at zero lag equals one. See \code{\link{computeEffectiveAutoCorr}}
#' @param effCov alternative to specifying effCor: numeric vector of 
#'  effective covariance components
#'  first entry is the variance. See \code{\link{computeEffectiveAutoCorr}}
#' @param nEff possibility to specify precomputed number of effective 
#'  observations for speedup.
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
  , effCor = if (missing(effCov)) computeEffectiveAutoCorr(x) else 
      effCov/var(x, na.rm = TRUE)
  , na.rm = FALSE 
  , effCov # compatibility to seCor up to version 0.1.8
  , nEff = computeEffectiveNumObs(x, effCor, na.rm = na.rm)
){
  n <- if (na.rm) length(na.omit(x)) else length(x)
  if (n < 2) return(NA_real_)
  if (var(x, na.rm = TRUE) == 0) return(0)
  varCorVal <- varCor(x, nEff = nEff, na.rm = TRUE)/nEff
  sqrt(varCorVal)
}

#' Compute the unbiased variance accounting for empirical autocorrelations
#'
#' @param x numeric vector
#' @param na.rm logical. Should missing values be removed?
#' @param effCor numeric vector of effective correlation components
#'  first entry at zero lag equals one. See \code{\link{computeEffectiveAutoCorr}}
#'  The effective correlation is passed to \code{\link{computeEffectiveNumObs}}.
#' @param nEff possibility to specify precomputed number of effective 
#'  observations for speedup.
#'  
#' @details The default uses empirical autocorrelation
#'  estimates from the supplied data up to first negative component.
#'  For short series of \code{x} it is strongly recommended to to
#'  provide \code{effCov} that was estimated on a longer time series.
#'
#' @export
#' @return numeric scalar of unbiased variation of x
varCor <- function(
  x  
  , effCor = computeEffectiveAutoCorr(x) 
  , na.rm = FALSE 
  , nEff = computeEffectiveNumObs(x, effAcf = effCor)
){
  var_uncorr = var(x, na.rm=na.rm) 
  if (!is.finite(var_uncorr) || (var_uncorr == 0)) return(var_uncorr)
  n = length(x)
  if (n <= 1) return(var_uncorr)
  # BLUE Var(x) for correlated: Zieba11 eq.(1) 
  nmiss = sum(is.na(x))
  nfin = n - nmiss
  var_uncorr*(nfin-1)*nEff/(nfin*(nEff-1))
}

seCor_depr2103 <- function(
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
  # note that we use autocovariance instead of autocorrelation
  varCorB <- (g0 + 2*sum( (n - k)/n * g))/n
  sqrt(varCorB)
}


#' Compute the effective number of observations taking into account autocorrelation
#' 
#' @param res numeric of autocorrelated numbers, usually observation -
#'  model residuals
#' @param effAcf autocorrelation coefficients.
#'  The first entry is fixed at 1 for zero distance.
#' @param na.rm if not set to TRUE will return NA in there are missings
#'  in the series
#' @param exact.na if set to FALSE then do not count and correct for missing
#'   in the sum of autocorrelation terms. This is faster, but results are
#'   increasingly biased high with increasing number of missings. 
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
#' @details Assumes records of all times present. 
#' DO NOT REMOVE OR FILTER NA records before. 
#' The length of the time series is used.
#' @details Handling of NA values: The formula from Zieba 2011 is extended
#' to subtract the number of missing pairs in the count of correlation terms.
#' If `exact.na=false` the original formula is used (after trimming edge-NAs).
#' 
#' @export
#' @return integer scalar: effective number of observations
#' @examples
#' # generate autocorrelated time series
#' res <- stats::filter(rnorm(1000), filter = rep(1,5), circular = TRUE)
#' res[100:120] <- NA
#' # plot the series of autocorrelated random variables
#' plot(res)
#' # plot their empirical autocorrelation function
#' acf(res, na.action = na.pass)
#' #effAcf <- computeEffectiveAutoCorr(res)
#' # the effective number of parameters is less than number of 1000 samples
#' (nEff <- computeEffectiveNumObs(res, na.rm = TRUE))
computeEffectiveNumObs <- function(
  res  
  , effAcf = computeEffectiveAutoCorr(res)
  , na.rm = FALSE 
  , exact.na = TRUE
){
  if (!isTRUE(na.rm) & any(is.na(res))) return(NA_integer_)
  resTr <- .trimNA(res)
  lacf = length(effAcf) -1 # acf starts with lag 0
  isFin <- is.finite(resTr)
  if (lacf < 1) return(sum(isFin))
  n <- length(resTr)
  if (n < 2) return(n)
  k = 1:min(n-1,lacf) 
  if (exact.na) {
    # number of missing combinations due to missing in x
    mka = count_NA_forlags(resTr, 0:length(k))
    m0 = mka[1]
    mk = mka[-1]
    nf = n - m0
    neff = nf/(1 + 2/nf*sum((n - k -mk) * effAcf[k+1]))  
  } else {
    neff = n/(1 + 2/n*sum((n - k) * effAcf[k+1]))  
  }
  neff
}

count_NA_forlags <- function(x, lags=0:(length(x)-1)) {
  lx = length(x)
  xb = is.na(x)
  cntna0 = sum(xb)
  vapply(lags, function(lag) {
    if (lag == 0) return(cntna0)
    if (lag > (lx-1)) stop(
      "expected lag < length(x)-1=",lx-1," but got lag ",lag)
    xm = cbind(xb[1:(lx-lag)],xb[-(1:lag)])
    sum(rowSums(xm) > 0)
  },0)
}

computeEffectiveNumObs_depr2103 <- function(
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
#' @examples
#' # generate autocorrelated time series
#' res <- stats::filter(rnorm(1000), filter = rep(1,5), circular = TRUE)
#' res[100:120] <- NA
#' (effAcf <- computeEffectiveAutoCorr(res))
computeEffectiveAutoCorr <- function(
  res, type = "correlation" ){
  # first compute empirical autocorrelations
  ans <- acf(res, na.action = na.pass, plot = FALSE, type = type)
  # next get the number of elements before crossing the zero line
  nC <- suppressWarnings(min(which(ans$acf <= 0)) - 1)
  if (!is.finite(nC)) {
    # if there was no below zero correlation within default lag.max then
    # repeat acf with computing all lags
    ans <- acf(res, na.action = na.pass, plot = FALSE, type = type, lag.max = Inf)
    # append -1 so that nC equals to full length if no negative correlation
    nC <- min(which(c(ans$acf,-1) <= 0)) - 1
  }
  nC <- pmax(1,nC) # return at least one component
  ans$acf[1:nC]
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

