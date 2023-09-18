#' Convert: \code{loom} ==> \code{SingleCellLoomExperiment}
#'
#' @inheritParams converters
#' @inheritDotParams LoomExperiment::SingleCellLoomExperiment
#' @export
#' @examples
#' library(Seurat)
#' obj <-  example_obj("loom")
#' sce <- loom_to_scle(obj)
loom_to_scle <- function(obj,
                         verbose=TRUE,
                         ...){
  requireNamespace("LoomExperiment")

  messager_to_()
  obj2 <- LoomExperiment::SingleCellLoomExperiment(obj,
                                                   ...)
  return(obj2)
}
