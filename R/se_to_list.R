#' Convert: \code{SummarizedExperiment} ==> \code{list}
#'
#' @export
#' @examples
#' obj <- example_obj("se")
#' obj2 <- se_to_list(obj)
se_to_list <- function(obj){
  list(data=obj@assays,
       obs=obj@colData,
       var=SummarizedExperiment::rowData(obj))
}
