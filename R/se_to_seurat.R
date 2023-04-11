#' Convert: \code{SummarizedExperiment} ==> \code{Seurat}
#'
#' @examples
#' obj <- example_obj("se")
#' obj2 <- se_to_seurat(obj)
se_to_seurat <- function(obj,
                         verbose=TRUE){
  messager("+ SummarizedExperiment ==> Seurat",v=verbose)
  obj <- se_to_sce(obj = obj,
                   verbose = verbose)
  Seurat::as.Seurat(obj)
}


