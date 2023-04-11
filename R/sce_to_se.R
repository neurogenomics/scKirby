#' Convert: \code{SingleCellExperiment} ==> \code{SummarizedExperiment}
#'
#' @export
#' @examples
#' obj <- example_obj("sce")
#' obj2 <- sce_to_se(obj)
sce_to_se <- function(obj,
                      verbose=TRUE){
  messager("+ SingleCellExperiment ==> SummarizedExperiment",v=verbose)
  methods::as(obj,"SummarizedExperiment")
}
