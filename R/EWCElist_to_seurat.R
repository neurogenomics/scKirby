

#' Convert: \code{EWCElist} ==> \code{Seurat}
#'
#' @examples
#' library(scKirby)
#' data("example_EWCElist")
#' seurat <- EWCElist_to_seurat(object=example_EWCElist)
EWCElist_to_seurat <- function(object,
                               verbose=T){
  messager("+ EWCElist ==> Seurat",v=verbose)
  seurat <- Seurat::CreateSeuratObject(counts = as(as.matrix(object$exp), "sparseMatrix"),
                                       meta.data = data.frame(object$annot, row.names = colnames(object$exp))
                                       )
  return(seurat)
}
