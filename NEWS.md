# lognorm 0.1.9
- Less bias with missing values in computeEffectiveNumObs
- New function varCor to compute unbiased variance of uncorrelated time series
- seCor now based on varCor, which reduces bias for small number of effective
  observations.

# lognorm 0.1.8
- Implement the Lo 2012 approximation of the distribution of the difference of
  two lognormally distributed random variables by a shifted lognormal 
  distribution.
- Add Formulas to vignette of lognormalSum  
- Add shift argument to moments, mode, and median to deal with shifted
  lognormal distribution.
- Implement sample-based Distribution function for the difference of two 
  lognormal variables.

# lognorm 0.1.7

- rewrite documentation using roxygen2 instead of inlinedocs
- getParmsLognormForLowerAndUpper: remove argument isTransScale and provide
  new function getParmsLognormForLowerAndUpperLog instead
- merge vignette aggregateCorrelated to vignette lognormalSum

# lognorm 0.1.6
fix for CRAN: reformat (html-escape) doi in DESCRIPTION

# lognorm 0.1.5

- Sum of lognormals: if all values are gap-filled, assume multiplicative
   standard deviation (sigmaStar) of the sum to be the mean of given sigmaStar.
   Before took it assumed to be the maximum. But this led to larger confidence 
   bounds compared to using the normal assumption.
   
- Estimate parameters from various statistics, such as mean and an
    upper quantile value, i.e. a practical maximum.

# lognorm 0.1.4

Support common case of estimating standard error in the presence of correlations:
function `seCor`.

# lognorm 0.1.3

estimateSumLognormal by default now returns NA if there are NA values in terms.
New argument na.rm = TRUE allows for previous behavior of neglecting those terms.

# lognorm 0.1.2

- Compute effective number of observations
  - fixed bracketing bug in formula
  - better dealing with NA at begin or end

# lognorm 0.1.1

basic setup
