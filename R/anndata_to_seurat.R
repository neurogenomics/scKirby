#' @describeIn anndata_to anndata_to
#'
#' @param ... Parameters passed to \code{sceasy::convertFormat}.
#' @export
#' @importFrom sceasy convertFormat
#' @examples
#' obj <- example_obj("anndata")
#' seurat <- anndata_to_seurat(obj)
anndata_to_seurat <- function(obj,
                              verbose=TRUE,
                              ...){
  messager_to()
  sceasy::convertFormat(obj = obj,
                        from = "anndata",
                        to = "seurat",
                        ...)
}
