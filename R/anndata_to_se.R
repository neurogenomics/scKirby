#' Convert: \code{Anndata} ==> \code{SummarizedExperiment}
#'
#' @inheritParams converters
#' @inheritParams to_se
#' @export
#' @examples
#' obj <- example_obj("anndata")
#' sce <- anndata_to_se(obj)
anndata_to_se <- function(obj,
                           as_sce=FALSE,
                           verbose=TRUE){
  messager("+ AnnData ==> SummarizedExperiment",v=verbose)
  obj2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(raw = DelayedArray::DelayedArray(
      methods::as( Matrix::t(obj$X), "sparseMatrix"))
      ),
    colData = obj$obs,
    rowData = obj$var
  )
  if(isFALSE(as_sce)){
    obj2 <- sce_to_se(obj = obj2,
                      verbose = verbose)
  }
  obj2 <- check_se_rownames(obj2, verbose = verbose)
  return(obj2)
}
