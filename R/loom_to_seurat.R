#' Convert: \code{loom} ==> \code{Seurat}
#'
#' @export
#' @examples
#' library(Seurat)
#' obj <- example_obj("loom")
#' obj2 <- loom_to_seurat(obj)
loom_to_seurat <- function(obj,
                           verbose=TRUE){
  messager("+ loom ==> Seurat",v=verbose)
  obj2 <- Seurat::as.Seurat(obj)
  return(obj2)
}
