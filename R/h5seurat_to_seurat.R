#' Convert: \code{Seurat} ==> \code{H5Seurat}
#'
#' @inheritParams converters
#' @inheritDotParams SeuratObject::as.Seurat
#'
#' @export
#' @examples
#' obj <- example_obj("h5seurat")
#' seurat <- h5seurat_to_seurat(obj)
h5seurat_to_seurat <- function(obj,
                               verbose = TRUE,
                               ...){

  messager("+ Seurat ==> h5Seurat",v=verbose)
  Seurat::as.Seurat(obj,
                    verbose = verbose,
                    ...)
}
