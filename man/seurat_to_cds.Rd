% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat_to_cds.R
\name{seurat_to_cds}
\alias{seurat_to_cds}
\title{Convert: \code{Seurat} ==> \code{CellDataSet}}
\usage{
seurat_to_cds(obj, verbose = TRUE)
}
\arguments{
\item{obj}{A single-cell object supported by \pkg{scKirby}.
See \link[scKirby]{converters} for a table of all supported conversions.}

\item{verbose}{Print messages.}
}
\description{
Convert: \code{Seurat} ==> \code{CellDataSet}
}
\examples{
obj <- example_obj("seurat")
obj2 <- seurat_to_cds(obj)
}
