% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{buildSignatureMatrixMAST}
\alias{buildSignatureMatrixMAST}
\title{#' Building the signature matrix using mast}
\usage{
buildSignatureMatrixMAST(
  scdata,
  id,
  path = NULL,
  verbose = FALSE,
  ncores = 1,
  diff.cutoff = 0.5,
  pval.cutoff = 0.01
)
}
\arguments{
\item{scdata}{The single cell data matrix}

\item{id}{A Vector of the cell type annotations}

\item{path}{OPTIONAL path for saving generated files}

\item{verbose}{Whether to produce an output on the console.}

\item{ncores}{How many cores to use for DGE analysis; NULL to use all available cores-1}

\item{diff.cutoff}{The FC cutoff}

\item{pval.cutoff}{The pValue cutoff}
}
\value{
The computed signature matrix
}
\description{
When path = NULL, the generated files in the processes will not be saved and output.
}
