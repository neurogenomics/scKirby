% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_anndata.R
\name{to_anndata}
\alias{to_anndata}
\title{To \code{AnnData}}
\usage{
to_anndata(
  obj,
  save_path = tempfile(fileext = ".h5ad"),
  reimport = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{obj}{A single-cell object supported by \pkg{scKirby}.
See \link[scKirby]{converters} for a table of all supported conversions.}

\item{save_path}{Path to save the \code{H5Seurat} to.}

\item{reimport}{Save and re-import the \code{AnnData} object into R to ensure
all data has been converted from Python-native to R-native objects
(e.g. pandas data.frames vs. R data.frames).}

\item{verbose}{Print messages.}
}
\description{
Convert any single-cell object to \code{AnnData} format.
}
\examples{
obj <- example_obj("cds")
obj2 <- to_anndata(obj)
}
