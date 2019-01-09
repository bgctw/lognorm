# lognorm 0.1.5

- Sum of lognormals: if all values are gap-filled, assume multiplicative
   standard deviation (sigmaStar) of the sum to be the mean of given sigmaStar.
   Before took it assumed to be the maximum. But this led to larger confidence 
   boundsas compared to using the normal assumption.

# lognorm 0.1.4

Support common case of estimating standard error in the presence of correlations:
function `seCor`.

# lognorm 0.1.3

EstimateSumLognormal by default now returns NA if there are NA values in terms.
New argument na.rm = TRUE allows for previous behaviour of neglecting those terms.

# lognorm 0.1.2

- Computate of effective number of observations
  - fixed bracketing bug in formula
  - better dealing with NA at begin or end

# lognorm 0.1.1

basic setup
