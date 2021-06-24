
#' Convert: \code{SingleCellExperiment} ==> \code{Seurat}
#'
#' @examples
#' library(scKirby)
#' data("example_sce")
#' seurat <- sce_to_seurat(example_sce)
#' @examples
sce_to_seurat <- function(object,
                          verbose=T){
  messager("+ SingleCellExperiment ==> Seurat",v=verbose)
  seurat <- Seurat::as.Seurat(object)
  return(seurat)
}


