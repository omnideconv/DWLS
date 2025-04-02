\name{trimData}
\alias{trimData}
\title{Trims bulk data and single-cell signature matrix to contain the same genes}
\description{Trims bulk data and single-cell signature matrix to contain the same genes}
\usage{
    trimData(Signature, bulkData)
}
\arguments{
    \item{Signature}{Gene by cell-type matrix consisting of averaged expression values across cells from each cell type}
    \item{bulkData}{List containing bulk data values.}
}
\value{
    List containing trimmed Signature (S) and bulk (B) matrices.
}
\examples{
    trimmed <- trimData(Signature,bulkData)
    S <- trimmed$S
    B <- trimmed$B
}
\keyword{geneFilter}
