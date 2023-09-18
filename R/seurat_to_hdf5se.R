#' Convert: \code{Seurat} ==> \code{hdf5se}
#'
#' @inheritDotParams save_hdf5se
#' @inheritParams converters
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- seurat_to_hdf5se(obj)
seurat_to_hdf5se <- function(obj,
                             save_path = tempfile(),
                             verbose = TRUE,
                             ...){

  messager("+ Seurat ==> hdf5se",v=verbose)
  obj <- seurat_to_se(obj = obj,
                      verbose = verbose)
  save_hdf5se(obj = obj,
              save_dir = save_path,
              verbose = verbose,
              ...)
}
