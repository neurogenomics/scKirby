

#' Subsetted mouse ESC \code{SingleCellExperiment}  object
#'
#' A dataset containing 300 cells and 2026 genes from two batches of mouse ESC data.
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
