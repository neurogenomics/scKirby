#' Convert: \code{SummarizedExperiment} ==> \code{SingleCellLoomExperiment}
#'
#' @param as_scle Convert to \link[LoomExperiment]{LoomExperiment} format.
#' @export
#' @examples
#' obj <- example_obj("se")
#' obj2 <- se_to_le(obj)
se_to_le <- function(obj,
                     as_scle=FALSE,
                     verbose=TRUE){
  messager("+ SummarizedExperiment ==> SingleCellLoomExperiment",v=verbose)
  obj2 <- sceasy::convertFormat(obj = obj, from = "sce", to = "loom")
  if(isFALSE(as_scle)){
    obj2 <- scle_to_le(obj = obj2,
                       verbose = verbose)
  }
  return(obj2)
}
