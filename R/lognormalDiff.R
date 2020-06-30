#' Estimatd the shifted-lognormal approximation to difference of two lognormals
#' 
#' The distribtuion of y = a - b + s, where a and b are two lognormal random
#' variables and s is a constant to be estimated, can be approximated
#' by a lognormal distribution.
#' 
#' @param mu_a center parameter of the first term
#' @param mu_b center parameter of the second term
#' @param sigma_a scale parameter of the first term
#' @param sigma_b scale parameter of the second term
#' @param corr correlation between the two random variables
#'
#' @return numeric vector with components mu, sigma, and shift, the components
#' of the shifted lognormal distribution.
#' @export
estimateDiffLognormal <- function(mu_a, mu_b, sigma_a, sigma_b, corr = 0){
  if (!(sigma_a > sigma_b)) stop(
    "expected sigma_a > sigma_b but got ", sigma_a, " <= ", sigma_b,
    ". Exchange the terms and negate the resulting quantiles ",
    "(see vignette('lognormalDiff').")
  Sa = exp(mu_a + sigma_a^2/2)
  Sb = exp(mu_b + sigma_b^2/2)
  sigma2_m <- sigma_a^2 + sigma_b^2 - 2*corr*sigma_a*sigma_b
  sigma_mt <- (sigma_a^2 - sigma_b^2)/(2*sqrt(sigma2_m))
  S0p = Sa + Sb #eq. 2.2.
  S0m = Sa - Sb
  if (S0m/S0p > 1e-2) warning(
    "Expected S0+/S0- << 1 but this ratio was ",S0m/S0p, ". The Lo 2012 ",
    "approximation becomes inaccurate for small number a and b.")
  shift = sigma2_m/(sigma_a^2 - sigma_b^2)*S0p
  S0mt <- S0m + shift
  if (sigma_mt > 0.2) warning(
    "Expected small sigma of the shifted distribution but got ",sigma_mt,
    ". The Lo 2012 approximation becomes inaccurate for larger sigma.")
  muSum <- log(S0mt - sigma_mt^2/2)
  return(c(mu = muSum, sigma = sigma_mt, shift = shift))
}

