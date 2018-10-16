\name{getParmsLognormForExpval}
\alias{getParmsLognormForExpval}
\title{getParmsLognormForExpval}
\description{get the lognormal parameters by expected value}
\usage{getParmsLognormForExpval(mean, sigmaStar)}
\arguments{
  \item{mean}{expected value at original scale}
  \item{sigmaStar}{multiplicative standard deviation}
}



\author{Thomas Wutzler}




\examples{
.mean <- 1
.sigmaStar <- c(1.3,2)
(parms <- getParmsLognormForExpval(.mean, .sigmaStar))
# multiplicative standard deviation must equal the specified value
cbind(exp(parms[,"sigma"]), .sigmaStar)
}
