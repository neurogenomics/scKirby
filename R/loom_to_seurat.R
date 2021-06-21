

#' Convert: \code{loom} ==> \code{Seurat}
#'
#' @examples
#' library(scKirby)
#' loom <- example_loom()
#' sce <- loom_to_seurat(loom)
#' @examples
loom_to_seurat <- function(object,
                           verbose=T){
  messager("+ loom ==> Seurat",v=verbose)
  seurat <- Seurat::as.Seurat(loom)
  return(seurat)
}
