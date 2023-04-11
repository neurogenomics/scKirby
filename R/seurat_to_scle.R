#' Convert: \code{Seurat} ==> \code{LoomExperiment}
#'
#' @inheritParams se_to_le
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- seurat_to_le(obj)
seurat_to_le <- function(obj,
                         as_scle=FALSE,
                         verbose=TRUE){
  obj2 <- seurat_to_se(obj = obj,
                       as_sce = TRUE,
                       verbose = verbose)
  obj2 <- se_to_le(obj = obj2,
                   as_scle = as_scle,
                   verbose = verbose)
  return(obj2)
}
