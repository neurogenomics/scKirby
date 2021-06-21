
#' Convert: \code{Anndata} ==> \code{SingleCellExperiment}
#'
#' @examples
#' library(scKirby)
#' data("example_anndata")
#' sce <- anndata_to_sce(example_anndata)
#' @examples
anndata_to_sce <- function(object,
                           verbose=T){
  messager("+ AnnData ==> SingleCellExperiment",v=verbose)
  # sceasy::convertFormat(obj = object, from = "anndata", to="seurat",
  #                       outFile = "~/Desktop/tmp.rds")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays      = list(raw = DelayedArray::DelayedArray(as( SparseM::t(object$X), "sparseMatrix"))),
    colData     = object$obs,
    rowData     = object$var
  )
  sce <- check_sce_rownames(sce, verbose = verbose)
  return(sce)
}
