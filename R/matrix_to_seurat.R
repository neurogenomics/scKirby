


#' Convert: \code{matrix} ==> \code{Seurat}
#'
#' @examples
#' library(scKirby)
#' data("example_seurat")
#' sce <- matrix_to_seurat(object=example_seurat@assays$RNA@counts)
matrix_to_seurat <- function(object,
                          verbose=T){
  messager("+ Matrix ==> Seurat",v=verbose)
  seurat <- Seurat::CreateSeuratObject(counts = object,
                                       meta.data = data.frame(cellid=colnames(object),
                                                              row.names = colnames(object)))
  return(seurat)
}
