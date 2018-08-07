#' @export
estimateSumLognormalSample <- function(
  ### Estimate the parameters of the lognormal approximation to the sum
  mu        ##<< numeric vector of center parameters of terms at log scale
  , sigma   ##<< numeric vector of variance parameter of terms at log scale
  , resLog  ##<< time series of model-residuals at log scale 
  ## to estimate correlation
  , effAcf = computeEffectiveAutoCorr(resLog) ##<< effective autocorrelation
  ## coefficients (may provide precomputed for efficiency or if the sample
  ## of resLog is too small) set to 1 to assume uncorrelated sample
  , isGapFilled = logical(0) ##<< logical vector whether entry is gap-filled 
  ## rather than an original measurement, see details
){
  # only one term, return the parameters
  if (length(mu) == 1 ) return(
    c(mu = as.vector(mu), sigma = as.vector(sigma), nEff = 1)
  )
  nEff <- computeEffectiveNumObs(resLog, effAcf = effAcf)
  ##details<<
  ## If there are no gap-filled values, i.e. \code{all(!isGapFilled)} or
  ## \code{!length(isGapFilled)} (the default), distribution parameters
  ## are estimated using all the samples. Otherwise, the scale parameter
  ## (uncertainty) is first estimated using only the non-gapfilled records.
  if (length(isGapFilled) && any(isGapFilled)) {
    isMeasured <- !isGapFilled
    sigmaSum <- if (!sum(isMeasured)) {
      ##details<< If there are only gap-filled records, assume uncertainty
      ## to be the largest uncertainty of given gap-filled records.
      max(sigma)
    } else {
      # recursive call with only the measured records
      ans1 <- estimateSumLognormalSample(
        mu[isMeasured], sigma[isMeasured], resLog[isMeasured]
        ,effAcf = effAcf, isGapFilled = logical(0)
      )
      ans1["sigma"]
    }
    p <- estimateSumLognormal(
      # sigma will not be used only checked for is.finite()
      mu, sigma = 0, sigmaSum = sigmaSum, corrLength = length(effAcf) - 1)
    return(c(p , nEff = nEff))
  }
  ##details<<
  ## Correlation matrix is constructed by inspecting the residuals of
  ## model predictions - observations at log-scale. By default, only the first auto
  ## correlation components are used that are positive
  ## (see \code{\link{computeEffectiveAutoCorr}}).
  corrM <- diag(nrow = length(mu))
  if (length(effAcf) > 1) {
    corrM <- setMatrixOffDiagonals(
      corrM, value = effAcf[-1], isSymmetric = TRUE)
  }
  ##value<< numeric vector with components "mu", "sigma", and "nEff"
  ## the parameters of the lognormal distribution at log scale
  ## (Result of \code{link{estimateSumLognormal}})
  ## and the number of effective observations.
  p <- estimateSumLognormal( mu, sigma, corr = corrM, corrLength = length(effAcf) - 1)
  return(c(p , nEff = nEff))
}


#' @export
estimateSumLognormal <- function(
  ### Estimate the distribution parameters of the lognormal approximation to the sum
  mu       ##<< numeric vector of center parameters of terms at log scale
  , sigma  ##<< numeric vector of variance parameter of terms at log scale
  , corr = diag(nrow = length(mu)) ##<< numeric matrix of correlations between
  ## the random variables
  , sigmaSum = numeric(0) ##<< numeric scalar: possibility to specify
  ## of a precomputed scale parameter
  , corrLength = nTerm   ##<< integer scalar: set correlation length to
  ## smaller values
  ## to speed up computation by neglecting correlations among terms
  ## further apart.
  ## Set to zero to omit correlations.
  , isStopOnNoTerm = FALSE ##<< if no finite estiamte is provide, by
  ## default, NA is returned for the sum.
  ## Set this to TRUE to issue an error instead.
  , effAcf                 ##<< numeric vector of effective autocorrelation
  ## This overides arguments \code{corr} and \code{corrLength}
){
  ##references<< 
  ## Lo C (2013) WKB approximation for the sum of two correlated lognormal 
  ## random variables.
  ## Applied Mathematical Sciences, Hikari, Ltd., 7 , 6355-6367 
  ## 10.12988/ams.2013.39511
  iFinite <- which( is.finite(mu) & is.finite(sigma))
  muFin <- mu[iFinite]
  sigmaFin <- sigma[iFinite]
  nTerm = length(muFin)
  if (nTerm == 0) if (isTRUE(isStopOnNoTerm)) {
    stop( "need to provide at least one term"
          ,", but length of finite mu and sigma was zero")
  } else {
    # if no finite term was given, return NA
    return(c(mu = NA_real_, sigma = NA_real_))
  }
  if (nTerm == 1) return(structure(c(muFin,sigmaFin), names = c("mu","sigma")))
  if (!missing(effAcf) && length(effAcf)) {
    corr <- getCorrMatFromAcf(length(mu), effAcf)
    corrLength <- length(effAcf)
  }
  corr <- corr[iFinite,iFinite]
  S = exp(muFin)
  Ssum = sum(S)
  sigma2Eff <- if (length(sigmaSum)) {
    # if sigma of the sum has been pre-specified
    sigmaSum^2
  } else if (corrLength == 0) {
      sumDiag <- sum(sigmaFin*sigmaFin*S*S)/Ssum^2
  } else {
    jStarts <- pmax(1, (1:nTerm) - corrLength)
    jEnds <- pmin(nTerm, (1:nTerm) + corrLength)
    ansi = sapply(1:nTerm, function(i){
      j <- jStarts[i]:jEnds[i]
      ansj <- corr[i,j]*sigmaFin[i]*sigmaFin[j]*S[i]*S[j]  #Ssum outside loop
      # ansj2 <- sapply(jStarts[i]:jEnds[i], function(j){
      #   #corr[i,j]*sigmaFin[i]*sigmaFin[j]*S[i]/Ssum*S[j]/Ssum
      #   corr[i,j]*sigmaFin[i]*sigmaFin[j]*S[i]*S[j]  #Ssum outside loop
      # })
      # sum(ansj2)
      sum(ansj)
    })
    sum(ansi)/Ssum^2
  }
  muSum = log(Ssum) + sigma2Eff/2
  ##value<< numeric vector with two components mu and sigma
  ## the parameters of the lognormal distribution at log scale
  return(c(mu = as.vector(muSum), sigma = as.vector(sqrt(sigma2Eff))))
}

estimateSumLognormalBenchmark <- function(
  ### Estimate the distribution parameters of the lognormal approximation to the sum
  mu      ##<< numeric vector of center parameters
  , sigma  ##<< numeric vector of variance parameter at log sclae
  , corr = diag(nrow = nTerm) # TODO correlations
){
  ##details<< Implements estimation according to
  ## Messica A(2016) A simple low-computation-intensity model for approximating
  ## the distribution function of a sum of non-identical lognormals for
  ## financial applications. 10.1063/1.4964963
  nTerm = length(mu)
  S = exp(mu)
  Ssum = sum(S)
  # diag
  ansi = sapply(1:nTerm, function(i){
    ansj = sapply(i, function(j){
      #corr[i,j]*sigma[i]*sigma[j]*S[i]/Ssum*S[j]/Ssum
      corr[i,j]*sigma[i]*sigma[j]*S[i]*S[j]
    })
    sum(ansj)
  })
  sumDiag = sum(ansi)
  #
  ansi = sapply(1:nTerm, function(i){
    ansj = sapply(1:nTerm, function(j){
      #corr[i,j]*sigma[i]*sigma[j]*S[i]/Ssum*S[j]/Ssum
      corr[i,j]*sigma[i]*sigma[j]*S[i]*S[j]
    })
    sum(ansj)
  })
  sigma2Eff <- sum(ansi)/Ssum^2
  # for performance reasons take out Ssum of the sum
  #sigma2Eff <- sum(ansi)/(nTerm*Ssum)^2
  muSum = log(Ssum) + sigma2Eff/2
  ##value<< numeric vector with two components mu and sigma
  ## the parameters of the lognormal distribution at log scale
  return(c(mu = muSum, sigma = sqrt(sigma2Eff)))
}



