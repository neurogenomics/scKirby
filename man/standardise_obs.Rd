% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardise_obs.R
\name{standardise_obs}
\alias{standardise_obs}
\title{Prepare metadata}
\usage{
standardise_obs(
  obj,
  save_path = NULL,
  species,
  standardise_species = TRUE,
  method = "gprofiler",
  level_cols = list(),
  extra_metadata = NULL,
  return_obs = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{obj}{A single-cell object supported by \pkg{scKirby}.
See \link[scKirby]{converters} for a table of all supported conversions.}

\item{save_path}{Path to save the modified \code{obj} to.}

\item{species}{The name of the species that the single-cell \code{obj}
came from. Can be provided as the scientific name (e.g. "Mus musculus") or
common name (e.g. "mouse").}

\item{standardise_species}{Standardise the name using
\link[orthogene]{map_species}.}

\item{method}{R package to use for gene mapping:
\itemize{
\item{\code{"gprofiler"} : Slower but more species and genes.}
\item{\code{"homologene"} : Faster but fewer species and genes.}
\item{\code{"babelgene"} : Faster but fewer species and genes.
Also gives consensus scores for each gene mapping based on a
 several different data sources.}
}}

\item{level_cols}{Names of metadata columns that indicate the hierarchically
organized cell type annotation levels.
Each level can be a combination of other columns, e.g.:
 \code{
   list(level1="organ",
        level2=c("celltype"),
        level3=c("organ","celltype"))
   }}

\item{extra_metadata}{Extra metadata to add to obs.
 Can be a:
 \itemize{
 \item{NULL : }{No extra metadata will be added (default).}
 \item{data.frame : }{With the same number of rows as observations (samples)
 in \code{obj}. NOTE: We assumed the rows of the \code{extra_metadata}
  data.frame are in the same order as the observations (samples)
 in \code{obj}}
 \item{list : }{A named list with a single entry per element,
 to be applied to all observation.}
}}

\item{return_obs}{Return just standardised \code{obs}
instead of the entire \code{obj}.}

\item{verbose}{Print messages.}
}
\description{
Standardise metadata for any single-cell object.
}
\examples{
obj <- example_obj()
obj2 <- standardise_obs()
}
