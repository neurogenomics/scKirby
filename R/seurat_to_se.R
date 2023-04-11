#' Convert: \code{Seurat} ==> \code{SummarizedExperiment}
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- seurat_to_se(obj)
seurat_to_se <- function(obj,
                          as_sce=TRUE,
                          verbose=TRUE){
  messager("+ Seurat ==> SingleCellExperiment",v=verbose)
  obj2 <- Seurat::as.SingleCellExperiment(obj)
  if(isFALSE(as_sce)){
    obj2 <- sce_to_se(obj = obj2,
                      verbose = verbose)
  }
  return(obj2)
}
