#' Convert: \code{Seurat} ==> \code{CellDataSet}
#'
#' @export
#' @examples
#' obj <- example_obj("Seurat")
#' obj2 <- seurat_to_cds(obj)
seurat_to_cds <- function(obj,
                          verbose=TRUE){
  Seurat::as.CellDataSet(obj)
}
