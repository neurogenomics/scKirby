#' Convert: \code{Seurat} ==> \code{H5Seurat}
#'
#' @param save_path Path to save the \code{H5Seurat} to.
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' seurat <- seurat_to_h5seurat(obj)
seurat_to_h5seurat <- function(obj,
                               save_path = tempfile(fileext = ".h5seurat"),
                               verbose = TRUE,
                               ...){

  messager("+ Seurat ==> h5Seurat",v=verbose)
  dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
  SeuratDisk::as.h5Seurat(obj,
                          filename = save_path,
                          overwrite = FALSE,
                          verbose = verbose,
                          ...)
}
