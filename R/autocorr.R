
#' @export
computeEffectiveNumObs <- function(
  ### compute the effective number of observations taking into account autocorrelation
  res  ##<< numeric of autocorrelated numbers, usually observation - model residuals
  , effAcf = computeEffectiveAutoCorr(res) ##<< may provide precomputed for efficiency
){
  n <- length(res)
  nC <- length(effAcf)
  if (!nC || (length(effAcf) == 1)) return(n)
  nEff <- n/(1 + 2*sum(1 - 1:nC/n*effAcf))
  ##value<< integer scalar: effective number of observations
  nEff
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
  ## Zieba 2011 Standard Deviation of the Mean of Autocorrelated
  ## Observations Estimated with the Use of the Autocorrelation Function Estimated
  ## From the Data
  n <- length(res)
  if (!is.finite(nC)) return(integer(0))
  ##value<< numeric vector: stongest compponents of the autocorrelation function
  ans$acf[1:nC]
}

#' @export
varEffective <- function(
  ### estimate the variance of a correlated time series
  res  ##<< numeric of autocorrelated numbers, usually observation - model residuals
  , nEff = computeEffectiveNumObs(res) ##<< effective number of observations
  , ...  ##<< further arguments to \code{\link{var}}
) {
  ##details<< The BLUE is not anymore the usual variance, but a modified
  ## variance as given in Zieba 2011
  a <- 1
  var(res, ...) * nEff/(nEff - 1)
  ### The estimated variance of the sample
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
){
  setMatrixOffDiagonals(
    diag(nrow = nRow), value = effAcf, isSymmetric = TRUE)
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
