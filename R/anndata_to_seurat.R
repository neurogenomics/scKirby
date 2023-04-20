#' Convert: \code{Anndata} ==> \code{Seurat}
#'
#' @export
#' @examples
#' obj <- example_obj("anndata")
#' seurat <- anndata_to_seurat(obj)
anndata_to_seurat <- function(obj,
                              verbose=TRUE,
                              ...){
  messager("+ AnnData ==> Seurat",v=verbose)
  sceasy::convertFormat(obj = obj,
                        from = "anndata",
                        to = "seurat",
                        ...)
}
