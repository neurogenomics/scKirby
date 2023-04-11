#' Convert: \code{SummarizedExperiment} ==> \code{SingleCellExperiment}
#'
#' @export
#' @examples
#' obj <- example_obj("se")
#' obj2 <- se_to_sce(obj)
se_to_sce <- function(obj,
                      verbose=TRUE){
  if(methods::is(obj,"SingleCellExperiment")) {
    obj2 <- check_se_rownames(obj,
                              verbose = verbose)
    return(obj2)
  }
  messager("+ SummarizedExperiment ==> SingleCellExperiment",v=verbose)
  obj2 <- methods::as(obj,"SingleCellExperiment")
  obj2 <- check_se_rownames(obj2, verbose = verbose)
  return(obj2)
}
