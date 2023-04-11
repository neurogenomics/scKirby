#' Convert: \code{Seurat} ==> \code{anndata}
#'
#' @param reimport Save and re-import the \code{anndata} object into R to ensure
#' all data has been converted from Python-native to R-native objects
#' (e.g. pandas data.frames vs. R data.frames).
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- seurat_to_anndata(obj)
seurat_to_anndata <- function(obj,
                              reimport=TRUE,
                              save_path= tempfile(fileext = ".h5ad"),
                              verbose=TRUE){
  adat <- sceasy::convertFormat(obj = obj,
                                from = "seurat",
                                to = "anndata")
  if(isTRUE(reimport)){
    adat$write_h5ad(filename = save_path)
    adat <- anndata::read_h5ad(save_path)
  }
  return(adat)
}
