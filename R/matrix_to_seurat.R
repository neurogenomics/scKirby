#' Convert: \code{matrix} ==> \code{Seurat}
#'
#' @inheritParams converters
#' @export
#' @examples
#' obj <- example_obj("matrix")
#' obj2 <- matrix_to_seurat(obj)
matrix_to_seurat <- function(obj,
                             verbose=TRUE){
  messager_to()
  obj2 <- Seurat::CreateSeuratObject(
    counts = obj,
    meta.data = data.frame(cellid = colnames(obj),
                           row.names = colnames(obj))
  )
  return(obj2)
}
