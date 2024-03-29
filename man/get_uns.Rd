% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_uns.R
\name{get_uns}
\alias{get_uns}
\title{Get unstructured metadata}
\usage{
get_uns(obj, n = NULL, verbose = TRUE)
}
\arguments{
\item{obj}{A single-cell object supported by \pkg{scKirby}.
See \link[scKirby]{converters} for a table of all supported conversions.}

\item{n}{The maximum number of elements to return from that list.}

\item{verbose}{Print messages.}
}
\value{
Named list of unstructured objects.
}
\description{
Extract unstructured metadata from any single-cell object.
}
\examples{
obj <- example_obj("seurat")
uns <- get_uns(obj)
}
