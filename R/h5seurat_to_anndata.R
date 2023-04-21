#' Convert: \code{h5Seurat} ==> \code{anndata}
#'
#' @inheritParams echoconda::activate_env
#' @inheritParams SeuratDisk::Convert
#' @inheritDotParams anndata::write_h5ad
#' @export
#' @examples
#' obj <- example_obj("h5seurat")
#' obj2 <- h5seurat_to_anndata(obj)
h5seurat_to_anndata <- function(obj,
                              reimport = TRUE,
                              save_path = tempfile(fileext = ".h5ad"),
                              conda_env = "r-reticulate",
                              overwrite = FALSE,
                              verbose=TRUE,
                              ...){
  SeuratDisk::Convert(obj$filename,
                      dest = save_path,
                      overwrite = overwrite,
                      verbose = verbose)
  adat <- read_anndata(path = save_path,
                       verbose = verbose,
                       conda_env = conda_env)
  #### Return ###
  return(adat)
}
