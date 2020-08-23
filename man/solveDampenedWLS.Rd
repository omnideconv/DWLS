\name{solveDampenedWLS}
\alias{solveDampenedWLS}
\title{Solve overdetermined system of equations, \eqn{Sx=B}, for x, using dampened weighted least squares.}
\description{
Dampened weighted least squares involves solving:
\deqn{min_{x>0} \sum_{i=1}^n w_i(B_i-(Sx)_i)^2}
through an iterative approach.
S is a gene by cell-type signature matrix, B is a list containing bulk data, and x is
the number of cells from each cell type present in the bulk data.
x is then normalized:
\deqn{\tilde{x} = \frac{x}{\sum_{j=1}^n x},} where j is the cell-type index.
}
\usage{
    solveDampenedWLS(S, B)   
}
\arguments{
    \item{S}{Gene by cell-type matrix consisting of averaged expression values across cells from each cell type }
    \item{B}{List containing bulk data values.}
}
\value{
    List containing proportion of each cell type present in bulk data set.
}
\examples{
    prop <- solveDampenedWLS(S,B)
}
\keyword{solveDampenedWLS}
