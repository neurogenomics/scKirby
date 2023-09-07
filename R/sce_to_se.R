#' Convert: \code{SingleCellExperiment} ==> \code{SummarizedExperiment}
#'
#' @inheritParams converters
#' @export
#' @examples
#' obj <- example_obj("sce")
#' obj2 <- sce_to_se(obj)
sce_to_se <- function(obj,
                      verbose=TRUE){
  messager("+ SingleCellExperiment ==> SummarizedExperiment",v=verbose)
  obj2 <- methods::as(obj,"SummarizedExperiment")
  #### Rownames don't get transferred correctly during conversion ####
  S4Vectors::rownames(obj2) <- S4Vectors::rownames(obj)
  return(obj2)
}
