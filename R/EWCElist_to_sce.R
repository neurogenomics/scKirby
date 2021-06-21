

#' Convert: \code{EWCElist} ==> \code{SingleCellExperiment}
#'
#' @examples
#' library(scKirby)
#' data("example_EWCElist")
#' sce <- EWCElist_to_sce(example_EWCElist)
EWCElist_to_sce <- function(object,
                            verbose=T){
  messager("+ EWCElist ==> SingleCellExperiment",v=verbose)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays      = list(raw = DelayedArray::DelayedArray(as(as.matrix(object$exp), "sparseMatrix"))),
    colData     =  object$annot,
    rowData     =  row.names(object$exp)
  )
  sce <- check_sce_rownames(sce, verbose = verbose)
  return(sce)
}
