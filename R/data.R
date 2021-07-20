

#' Example \code{SingleCellExperiment}
#'
#' A subsetted dataset containing 300 cells and 2026 genes from two batches of mouse ESC data.
#'
#' Copied from the R package \pkg{scMerge}.
#'
#' @format A \code{SingleCellExperiment} object
#' @source \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2600/}
#' @references Kolodziejczyk et al.
#' @examples
#' \dontrun{
#' if(!require(scMerge))BiocManager::install("scMerge")
#' data('example_sce', package = 'scMerge')
#' batch_names<-unique(example_sce$batch)
#' usethis::use_data(example_sce, overwrite = T)
#' }
#' @export
"example_sce"



#' Example \code{Seurat}
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' data("pbmc_small")
#' example_seurat <- Seurat::UpdateSeuratObject(pbmc_small)
#' usethis::use_data(example_seurat, overwrite = T)
#' }
#' @export
"example_seurat"



#' Example \code{EWCElist}
#'
#' @examples
#' \dontrun{
#' example_EWCElist <- ewceData::cortex_mrna()
#' usethis::use_data(example_EWCElist, overwrite = T)
#' }
#' @export
"example_EWCElist"


#' Example \code{CellTypeDataset}
#'
#' @examples
#' \dontrun{
#' example_ctd <- ewceData::ctd()
#' ### Remove dendrograms to reduce file size
#' example_ctd[[1]]$plotting <- NULL
#' example_ctd[[2]]$plotting <- NULL
#' usethis::use_data(example_ctd, overwrite = T)
#' }
#' @export
"example_ctd"


#### Removing for now, bc it makes scKirby require monocle (which can be hard to install)
##' Example \code{CellDataSet}
##'
##' @examples
##' \dontrun{
##' library(Seurat)
##' data("pbmc_small")
##' example_seurat <- Seurat::UpdateSeuratObject(pbmc_small)
##' example_cds <- Seurat::as.CellDataSet(example_seurat)
##' usethis::use_data(example_cds, overwrite = T)
##' }
##' @export
#"example_cds"



