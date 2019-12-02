#' @details
#' Essential functions are 
#' 
#' \itemize{
#' \item Compute summary statistics: 
#'   \code{\link{getLognormMoments}}
#' \item computing distribution parameters from summary statistics: 
#'   \code{\link{getParmsLognormForModeAndUpper}}
#' \item estimate distribution parameters from summary sample: 
#'   \code{\link{estimateParmsLognormFromSample}}
#' \item Approximate the sum of correlated lognormals: 
#'   \code{\link{estimateSumLognormalSample}}
#' }
#' 
#' Utilities for correlated data. These functions maybe moved to
#' a separate package in future.
#' \itemize{
#' \item Estimate summary statistics of autocorrelated data
#' \itemize{
#'   \item standard error of the mean: \code{\link{seCor}} 
#'   \item effective number of observations \code{\link{computeEffectiveNumObs}} 
#'   \item variance: \code{\link{varEffective}} 
#'   }
#' \item Return the vector of effective components of the autocorrelation: 
#'   \code{\link{computeEffectiveAutoCorr}} 
#' }
#' 
#' Otherwise refer to the vignettes 
#' @keywords internal
"_PACKAGE"
#> [1] "_PACKAGE"
#