#' Convert: \code{loom} ==> \code{SingleCellLoomExperiment}
#'
#' @export
#' @examples
#' library(Seurat)
#' obj <-  example_obj("loom")
#' sce <- loom_to_scle(obj)
loom_to_scle <- function(obj,
                         verbose=TRUE,
                         ...){
  requireNamespace("LoomExperiment")

  messager("+ loom ==> SingleCellLoomExperiment",v=verbose)
  obj2 <- LoomExperiment::SingleCellLoomExperiment(obj)
  return(obj2)
}
