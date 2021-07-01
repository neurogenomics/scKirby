

#' Convert: \code{CellDataSet} ==> \code{SingleCellExperiment}
#'
#' @examples
#' library(scKirby)
#' data("example_cds")
#' sce <- cds_to_sce(object=example_cds)
#' @examples
cds_to_sce <- function(object,
                       verbose=T,
                       as_sparse=T,
                       as_DelayedArray=F){
  messager("+ CellDataSet ==> SingleCellExperiment",v=verbose)
  X <- Biobase::exprs(object)
  obs <- Biobase::pData(object)
  var <- Biobase::fData(object)
  if(as_sparse) X <- as(X,"sparseMatrix")
  if(as_DelayedArray) X <- DelayedArray::DelayedArray(X)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays      = list(raw = X),
    colData     = obs,
    rowData     = var
  )
  sce <- check_sce_rownames(sce, verbose = verbose)
  return(sce)
}
