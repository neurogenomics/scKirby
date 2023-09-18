#' Convert: \code{list} ==> \code{SummarizedExperiment}
#'
#' @inheritParams converters
#' @inheritParams to_se
#' @export
#' @examples
#' obj <- example_obj("list")
#' obj2 <- list_to_se(obj)
#' obj3 <- list_to_se(obj, as_sce=TRUE)
list_to_se <- function(obj,
                       as_sce=FALSE,
                       verbose=TRUE){

  #### Distinguish data with the same or different number of rows ####
  nrows <- sapply(obj$data, nrow)
  same_dim <- nrows==nrows[1]
  diff_dim <- nrows!=nrows[1]
  #### Construct objects ####
  if(isTRUE(as_sce)){
    messager("Converting: list ==> SingleCellExperiment",v=verbose)
    if(sum(diff_dim)>0){
      altExps <- lapply(obj$data[diff_dim],matrix_to_se, verbose=verbose)
    } else {
      altExps <- list()
    }
    obj2 <- SingleCellExperiment::SingleCellExperiment(
      assays = obj$data[same_dim],
      rowData = if(is.null(obj$var)) list() else obj$var,
      colData = if(is.null(obj$obs)) list() else obj$obs,
      altExps = altExps)
  } else {
    messager("Converting: list ==> SummarizedExperiment",v=verbose)
    obj2 <- SummarizedExperiment::SummarizedExperiment(
      assays = obj$data[same_dim],
      rowData = if(is.null(obj$var)) list() else obj$var,
      colData = if(is.null(obj$obs)) list() else obj$obs)
  }
  return(obj2)

}
