#' Convert: \code{SingleCellLoomExperiment} ==> \code{LoomExperiment}
#'
#' @inheritParams converters
#' @export
#' @examples
#' obj <- example_obj("se")
#' obj2 <- se_to_le(obj)
scle_to_le <- function(obj,
                       verbose=TRUE){

  messager("+ SingleCellLoomExperiment ==> LoomExperiment",v=verbose)
  obj2 <- methods::as(obj,"LoomExperiment")
  #### methods only transfers first assay ####
  assays <- SummarizedExperiment::assayNames(obj)
  if(length(assays)>length(SummarizedExperiment::assayNames(obj2))){
    for(i in seq_len(length(assays))[-1]){
      SummarizedExperiment::assay(obj2,i = i, withDimnames = FALSE) <-
        SummarizedExperiment::assay(obj,i = i)
    }
  }
  #### methods doesn't transfer names ####
  SummarizedExperiment::assayNames(obj2) <- assays
  #### Add rownames back ####
  S4Vectors::rownames(obj2) <- S4Vectors::rownames(obj)
  return(obj2)
}
