#' Convert: \code{SummarizedExperiment} ==> \code{SingleCellLoomExperiment}
#'
#' @export
#' @examples
#' obj <- example_obj("se")
#' obj2 <- se_to_scle(obj)
se_to_scle <- function(obj,
                       verbose=TRUE){
  messager("+ SummarizedExperiment ==> SingleCellLoomExperiment",v=verbose)
  sceasy::convertFormat(obj = obj, from = "sce", to = "loom")
}
