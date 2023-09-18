#' Convert: \code{Anndata} ==> \code{Seurat}
#'
#' NOTE: Do NOT try to \code{inheritDotParams sceasy::convertFormat}.
#' The documentation in that function is empty and this will a
#' cryptic error with \link[devtools]{document}.
#' @param ... Parameters passed to \code{sceasy::convertFormat}.
#' @inheritParams converters
#'
#' @export
#' @importFrom sceasy convertFormat
#' @examples
#' obj <- example_obj("anndata")
#' seurat <- anndata_to_seurat(obj)
anndata_to_seurat <- function(obj,
                              verbose=TRUE,
                              ...){
  messager_to_()
  sceasy::convertFormat(obj = obj,
                        from = "anndata",
                        to = "seurat",
                        ...)
}
