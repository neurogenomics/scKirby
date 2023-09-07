#' Convert: \code{Seurat} ==> \code{CellDataSet}
#'
#' @inheritParams converters
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- seurat_to_cds(obj)
seurat_to_cds <- function(obj,
                          verbose=TRUE){
  Seurat::as.CellDataSet(obj)
}
