#' Convert: \code{CellDataSet} ==> \code{Seurat}
#'
#' @inheritParams converters
#' @export
#' @examples
#' obj <- example_obj("cds")
#' obj2 <- cds_to_seurat(obj)
cds_to_seurat <- function(obj,
                          verbose=TRUE,
                          ...){
  messager_to()
  Seurat::as.Seurat(obj)
}
