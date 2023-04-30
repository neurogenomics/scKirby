#' Convert: \code{CellDataSet} ==> \code{AnnData}
#'
#' @inheritParams seurat_to_anndata
#' @export
#' @examples
#' obj <- example_obj("cds")
#' obj2 <- cds_to_anndata(obj)
cds_to_anndata <- function(obj,
                           reimport=TRUE,
                           save_path= tempfile(fileext = ".h5ad"),
                           verbose=TRUE){

  obj2 <- cds_to_seurat(obj = obj,
                        verbose = verbose)
  seurat_to_anndata(obj = obj2,
                    reimport = reimport,
                    verbose = verbose)
}
