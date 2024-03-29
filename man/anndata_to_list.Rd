% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anndata_to_list.R
\name{anndata_to_list}
\alias{anndata_to_list}
\title{Convert: \code{AnnData} ==> \code{list}}
\usage{
anndata_to_list(obj, verbose = TRUE)
}
\arguments{
\item{obj}{A single-cell object supported by \pkg{scKirby}.
See \link[scKirby]{converters} for a table of all supported conversions.}

\item{verbose}{Print messages.}
}
\description{
Convert: \code{AnnData} ==> \code{list}
}
\examples{
obj <- example_obj("anndata")
obj2 <- anndata_to_list(obj)
}
