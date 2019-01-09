\name{getParmsLognormForModeAndUpper}
\alias{getParmsLognormForModeAndUpper}
\title{getParmsLognormForModeAndUpper}
\description{Calculates mu and sigma of lognormal from mode and upper quantile.}
\usage{getParmsLognormForModeAndUpper(mle, quant, 
    sigmaFac = qnorm(0.99))}
\arguments{
  \item{mle}{numeric vector: mode at the original scale}
  \item{quant}{numeric vector: value at the upper quantile,
i.e. practical maximum}
  \item{sigmaFac}{sigmaFac=2 is 95\% sigmaFac=2.6 is 99\% interval}
}

\value{numeric matrix: columns mu and sigma parameter of the
lognormal distribution.
Rows correspond to rows of mle and quant}

\author{Thomas Wutzler}




\examples{
# example 1: a distribution with mode 1 and upper bound 5
(thetaEst <- getParmsLognormForModeAndUpper(1,5))
mle <- exp(thetaEst[1] - thetaEst[2]^2)
all.equal(mle, 1, check.attributes = FALSE)

# plot the distributions
xGrid = seq(0,8, length.out = 81)[-1]
dxEst <- dlnorm(xGrid, meanlog = thetaEst[1], sdlog = thetaEst[2])
plot( dxEst~xGrid, type = "l",xlab = "x",ylab = "density")
abline(v = c(1,5),col = "gray")

# example 2: true parameters, which should be rediscovered
theta0 <- c(mu = 1, sigma = 0.4)
mle <- exp(theta0[1] - theta0[2]^2)
perc <- 0.975		# some upper percentile, proxy for an upper bound
quant <- qlnorm(perc, meanlog = theta0[1], sdlog = theta0[2])
(thetaEst <- getParmsLognormForModeAndUpper(mle,quant = quant,sigmaFac = qnorm(perc)) )

#plot the true and the rediscovered distributions
xGrid = seq(0,10, length.out = 81)[-1]
dx <- dlnorm(xGrid, meanlog = theta0[1], sdlog = theta0[2])
dxEst <- dlnorm(xGrid, meanlog = thetaEst[1], sdlog = thetaEst[2])
plot( dx~xGrid, type = "l")
#plot( dx~xGrid, type = "n")
#overplots the original, coincide
lines( dxEst ~ xGrid, col = "red", lty = "dashed")

# example 3: explore varying the uncertainty (the upper quantile)
x <- seq(0.01,1.2,by = 0.01)
mle = 0.2
dx <- sapply(mle*2:8,function(q99){
  theta = getParmsLognormForModeAndUpper(mle,q99,qnorm(0.99))
  #dx <- dDistr(x,theta[,"mu"],theta[,"sigma"],trans = "lognorm")
  dx <- dlnorm(x,theta[,"mu"],theta[,"sigma"])
})
matplot(x,dx,type = "l")
}
