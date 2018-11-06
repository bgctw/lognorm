\name{seCor}
\alias{seCor}
\title{seCor}
\description{compute the standard error accounting for empirical correlations}
\usage{seCor(x)}
\arguments{
  \item{x}{numeric vector}
}
\details{computation according to 
https://stats.stackexchange.com/questions/274635/calculating-error-of-mean-of-time-series
but with considering only correlations up to first negative}
\value{numeric scalar of standard error of the mean of x}

\author{Thomas Wutzler}





