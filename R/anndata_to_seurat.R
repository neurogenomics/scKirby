#' Convert: \code{Anndata} ==> \code{Seurat}
#'
#'
#' @inheritParams converters
#' @param ... Parameters passed to \code{sceasy::convertFormat}.
#' @source \code{
#' NOTE: Do NOT try to "inheritDotParams sceasy::convertFormat".
#' The documentation in that function is empty and this will a
#' cryptic error with  devtools::document().
#' }
#'
#' @export
#' @importFrom sceasy convertFormat
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
