\name{DEAnalysisMAST}
\alias{DEAnalysisMAST}
\title{Calculates differentially expressed genes between each cell type and the rest for a given single-cell data set. (Using MAST function directly, without using Seurat)}
\description{
Finds differentially expressed genes between cell types using the MAST's lrTest function.
Saves differential expression information to specified path, in the form of a matrix 
for each cell type. This is used by the buildSignatureMatrixMAST function in DWLS.
}
\usage{
    DEAnalysisMAST(scdata, id, path)   
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
    DEAnalysisMAST(scdata, id, path)
}
\keyword{DEAnalysisMAST}
