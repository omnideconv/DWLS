\name{buildSignatureMatrixUsingSeurat}
\alias{buildSignatureMatrixUsingSeurat}
\title{Builds a signature matrix from single-cell data by calling Seurat.}
\description{
Given a single-cell matrix and cell type ids, the function constructs a gene by cell 
type matrix containing averaged expression values for each cell type, for highly
differentially expressed genes.  Differential expression analysis is conducted through Seurat.
}
\usage{
    buildSignatureMatrixUsingSeurat(scdata, id, path, diff.cutoff=0.5, pval.cutoff=0.01)
}
\arguments{
    \item{scdata}{Gene by cell single-cell data matrix.}
    \item{id}{List containing cell type ids corresponding to each column of the \code{scdata}}
    \item{path}{String corresponding to path in which to save DE information.}
    \item{diff.cutoff}{Scalar corresponding to the log fold change cutoff to use when picking differentially expressed genes (default=0.5)}
    \item{pval.cutoff}{Scalar corresponding to the pvalue cutoff to use when picking differentially expressed genes (default=0.01)}
}
\value{
    Gene by cell type matrix containing the gene expression signature of each cell type.
}
\examples{
    buildSignatureMatrixUsingSeurat(scdata, id, path, diff.cutoff=0.5, pval.cutoff=0.01)
}
\keyword{buildSignatureMatrixUsingSeurat}
