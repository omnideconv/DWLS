\name{DEAnalysis}
\alias{DEAnalysis}
\title{Calculates differentially expressed genes between each cell type and the rest for a given single-cell data set. (Using Seurat's function)}
\description{
Finds differentially expressed genes between cell types using the MAST test in Seurat.
Saves differential expression information to specified path, in the form of a matrix 
for each cell type. This is used by the buildSignatureMatrixUsingSeurat function in DWLS.
}
\usage{
    DEAnalysis(scdata, id, path)   
}
\arguments{
    \item{scdata}{Gene by cell single-cell data matrix.}
    \item{id}{List containing cell type ids corresponding to each column of the \code{scdata}}
    \item{path}{String corresponding to path in which to save DE information.}
}
\value{
    Saves DE information to specified path, in the form of a matrix for each cell type.
}
\examples{
    DEAnalysis(scdata, id, path)
}
\keyword{DEAnalysis}
