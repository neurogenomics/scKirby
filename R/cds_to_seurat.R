#' Convert: \code{CellDataSet} ==> \code{Seurat}
#'
#' @examples
#' obj <- example_obj("cds")
#' obj2 <- cds_to_seurat(obj)
cds_to_seurat <- function(obj,
                          verbose=TRUE){
  messager("+ CellDataSet ==> Seurat",v=verbose)
  Seurat::as.Seurat(obj)
}
