


#' Convert: \code{matrix} ==> \code{SingleCellExperiment}
#'
#' @examples
#' library(scKirby)
#' data("example_seurat")
#' sce <- matrix_to_sce(example_seurat@assays$RNA@counts)
matrix_to_sce <- function(object,
                          verbose=T){
  messager("+ Matrix ==> SingleCellExperiment",v=verbose)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays      = list(raw = DelayedArray::DelayedArray(as(as.matrix(object), "sparseMatrix"))),
  )
  sce <- check_sce_rownames(sce, verbose = verbose)
  return(sce)
}
