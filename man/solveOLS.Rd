\name{solveOLS}
\alias{solveOLS}
\title{Solve overdetermined system of equations, Sx=B, for x, using constrained ordinary least squares}
\description{
Constrained ordinary least squares is implemented using a quadratic programming 
approach, as such:
\deqn{min_{x>0}||B-Sx||^2}
S is a gene by cell-type signature matrix, B is a list containing bulk data, and x is
the number of cells from each cell type present in the bulk data.  x is then
normalized:
\deqn{\tilde{x} = \frac{x}{\sum_{j=1}^n x},} where j is the cell-type index.
}
\usage{
    solveOLS(S, B)   
}
\arguments{
    \item{S}{Gene by cell-type matrix consisting of averaged expression values across cells from each cell type }
    \item{B}{List containing bulk data values.}
}
\value{
    List containing proportion of each cell type present in bulk data set.
}
\examples{
    prop <- solveOLS(S,B)
}
\keyword{solveOLS}
