#' Convert: \code{CellDataSet} ==> \code{AnnData}
#'
#' @inheritParams seurat_to_anndata
#' @inheritParams echoconda::activate_env
#' @export
#' @examples
#' obj <- example_obj("cds")
#' obj2 <- cds_to_anndata(obj)
cds_to_anndata <- function(obj,
                           reimport=TRUE,
                           save_path= tempfile(fileext = ".h5ad"),
                           conda_env = "r-reticulate",
                           verbose=TRUE){

  obj2 <- cds_to_seurat(obj = obj,
                        verbose=verbose)
  seurat_to_anndata(obj = obj2,
                    reimport = reimport,
                    conda_env = conda_env,
                    verbose = verbose)
}
