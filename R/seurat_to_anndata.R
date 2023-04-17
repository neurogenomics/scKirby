#' Convert: \code{Seurat} ==> \code{anndata}
#'
#' @param reimport Save and re-import the \code{anndata} object into R to ensure
#' all data has been converted from Python-native to R-native objects
#' (e.g. pandas data.frames vs. R data.frames).
#' @inheritParams echoconda::activate_env
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- seurat_to_anndata(obj)
seurat_to_anndata <- function(obj,
                              reimport=TRUE,
                              save_path= tempfile(fileext = ".h5ad"),
                              conda_env = "r-reticulate",
                              verbose=TRUE){

  echoconda::activate_env(conda_env = conda_env,
                          method = "reticulate",
                          verbose = verbose)

  adat <- sceasy::convertFormat(obj = obj,
                                from = "seurat",
                                to = "anndata")
  if(isTRUE(reimport)){
    dir.create(dirname(save_path),showWarnings = FALSE, recursive = TRUE)
    # IMPORTANT! python cannot interpret "~"
    save_path <- path.expand(save_path)
    adat$write_h5ad(filename = save_path)
    adat <- anndata::read_h5ad(save_path)
  }
  return(adat)
}
