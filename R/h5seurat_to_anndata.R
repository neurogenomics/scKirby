#' Convert: \code{h5Seurat} ==> \code{anndata}
#'
#' @inheritParams converters
#' @inheritParams to_anndata
#' @inheritParams SeuratDisk::Convert
#' @inheritDotParams read_anndata
#' @export
#' @examples
#' obj <- example_obj("h5seurat")
#' obj2 <- h5seurat_to_anndata(obj)
h5seurat_to_anndata <- function(obj,
                                reimport = TRUE,
                                save_path = tempfile(fileext = ".h5ad"),
                                overwrite = FALSE,
                                assay=NULL,
                                verbose = TRUE,
                                ...){
  messager_to_()
  SeuratDisk::Convert(obj$filename,
                      dest = save_path,
                      overwrite = overwrite,
                      assay = assay,
                      verbose = verbose,
                      ...)
  adat <- read_anndata(path = save_path,
                       verbose = verbose,
                       ...)
  #### Return ###
  return(adat)
}
