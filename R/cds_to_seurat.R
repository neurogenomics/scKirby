

#' Convert: \code{CellDataSet} ==> \code{Seurat}
#'
#' @examples
#' library(scKirby)
#' data("example_cds")
#'
#' seurat <- cds_to_seurat(example_cds)
#' @examples
cds_to_seurat <- function(object,
                           verbose=T){
  messager("+ CellDataSet ==> Seurat",v=verbose)
  seurat <- Seurat::as.Seurat(object)
  return(seurat)
}
