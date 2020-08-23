\name{findDampeningConstant}
\alias{findDampeningConstant}
\title{Find a dampening constant for dampened weighted least squares using cross-validation.}
\description{
For a subset of genes, the function calculates the dampened weighted least squares 
solution for a particular dampening constant.  This is repeated 100 times per
dampening constant over a range of dampening constants.  The dampening constant with 
the lowest cross-validation variance is chosen.
}
\usage{
    findDampeningConstant(S, B, goldStandard)   
}
\arguments{
    \item{S}{Gene by cell-type matrix consisting of averaged expression values across cells from each cell type }
    \item{B}{List containing bulk data values.}
    \item{goldStandard}{An initial estimate of the number of each cell type present in the bulk data set. Often this estimate is derived from the \code{solveOLS} function.}
}
\value{
    Positive integer defining the dampening constant. The maximum weight value allowed will be \eqn{d=2^(j-1)}.
}
\examples{
    findDampeningConstant(S,B,solveOLS(S,B))
}
\keyword{findDampeningConstant}
