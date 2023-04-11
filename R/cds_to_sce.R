#' Convert: \code{CellDataSet} ==> \code{SummarizedExperiment}
#'
#' @export
#' @examples
#' obj <- example_obj("cds")
#' obj2 <- cds_to_se(obj)
cds_to_se <- function(obj,
                      as_sce=FALSE,
                      as_sparse=TRUE,
                      as_DelayedArray=FALSE,
                      verbose=TRUE){
  messager("+ CellDataSet ==> SummarizedExperiment",v=verbose)
  X <- Biobase::exprs(obj)
  obs <- Biobase::pData(obj)
  var <- Biobase::fData(obj)
  if(isTRUE(as_sparse)) X <- as(X,"sparseMatrix")
  if(isTRUE(as_DelayedArray)) X <- DelayedArray::DelayedArray(X)

  obj2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(raw = X),
    colData = obs,
    rowData = var
  )
  obj2 <- check_se_rownames(obj2, verbose = verbose)
  return(obj2)
}
