\name{setMatrixOffDiagonals}
\alias{setMatrixOffDiagonals}
\title{setMatrixOffDiagonals}
\description{set off-diagonal values of the matrix}
\usage{setMatrixOffDiagonals(x, diag = 1:length(value), 
    value, isSymmetric = FALSE)}
\arguments{
  \item{x}{numeric square matrix}
  \item{diag}{integer vector specifying the diagonals
0 is the center +1 the first
row to upper and -2 the second row to lower}
  \item{value}{numeric vector of values to fill in}
  \item{isSymmetric}{set to TRUE to to only
specify the upper diagonal element but also
set the lower in the mirrored diagonal}
}

\value{matrix with modified diagonal elements}

\author{Thomas Wutzler}





