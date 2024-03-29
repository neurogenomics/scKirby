% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/se_to_anndata.R
\name{se_to_anndata}
\alias{se_to_anndata}
\title{Convert: \code{SummarizedExperiment} ==> \code{AnnData}}
\usage{
se_to_anndata(
  obj,
  method = c("zellkonverter", "sceasy"),
  reimport = TRUE,
  save_path = tempfile(fileext = ".h5ad"),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{obj}{A single-cell object supported by \pkg{scKirby}.
See \link[scKirby]{converters} for a table of all supported conversions.}

\item{method}{Which R package to use for the conversion.}

\item{reimport}{Save and re-import the \code{AnnData} object into R to ensure
all data has been converted from Python-native to R-native objects
(e.g. pandas data.frames vs. R data.frames).}

\item{save_path}{Path to save the \code{H5Seurat} to.}

\item{verbose}{Print messages.}

\item{...}{Additional arguments passed to the respective converter function.}
}
\value{
\link[anndata]{AnnData} object.
}
\description{
Convert: \code{SummarizedExperiment} ==> \code{AnnData}
}
\examples{
obj <- example_obj("se")
obj2 <- se_to_anndata(obj)
}
