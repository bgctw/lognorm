\name{getLognormMode}
\alias{getLognormMode}
\title{getLognormMode}
\description{get the mode of a log-normal distribution}
\usage{getLognormMode(mu, sigma)}
\arguments{
  \item{mu}{center parameter (mean at log scale, log(median))}
  \item{sigma}{scale parameter (standard deviation at log scale)}
}

\value{the mode}

\author{Thomas Wutzler}




\examples{
# with larger sigma, the distribution is more skewed
# with mode further away from median = 1
getLognormMode(mu = log(1), sigma = c(log(1.2),log(2)))
}
