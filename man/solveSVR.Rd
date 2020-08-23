\name{solveSVR}
\alias{solveSVR}
\title{Solve overdetermined system of equations, \eqn{Sx=B}, for x, using nu-support
vector regression with a linear kernel}
\description{
Solve overdetermined system of equations, \eqn{Sx=B}, for x, using nu-support vector regression with a linear kernel.
S is a gene by cell-type signature matrix, B is a list containing bulk data, and x is
the number of cells from each cell type present in the bulk data.  x is then
normalized:
\deqn{\tilde{x} = \frac{x}{\sum_{j=1}^n x},} where j is the cell-type index.
}
\usage{
    solveSVR(S, B)   
}
\arguments{
    \item{S}{Gene by cell-type matrix consisting of averaged expression values across cells from each cell type }
    \item{B}{List containing bulk data values.}
}
\value{
    List containing proportion of each cell type present in bulk data set.
}
\examples{
    prop <- solveSVR(S,B)
}
\keyword{solveSVR}
