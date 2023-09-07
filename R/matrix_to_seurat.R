#' Convert: \code{matrix} ==> \code{Seurat}
#'
#' @inheritParams converters
#' @export
#' @examples
#' obj <- example_obj("matrix")
#' obj2 <- matrix_to_seurat(obj)
matrix_to_seurat <- function(obj,
                          verbose=TRUE){
  messager("+ Matrix ==> Seurat",v=verbose)
  obj2 <- Seurat::CreateSeuratObject(
    counts = obj,
    meta.data = data.frame(cellid=colnames(object),
                           row.names = colnames(object))
  )
  return(obj2)
}
