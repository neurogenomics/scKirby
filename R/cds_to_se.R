#' Convert: \code{CellDataSet} ==> \code{SummarizedExperiment}
#'
#' @inheritParams converters
#' @inheritParams to_se
#' @inheritParams get_x
#' @export
#' @examples
#' obj <- example_obj("cds")
#' obj2 <- cds_to_se(obj)
cds_to_se <- function(obj,
                      as_sce=FALSE,
                      as_sparse=TRUE,
                      as_delayedarray=FALSE,
                      verbose=TRUE){

  messager_to_()
  obj1 <- cds_to_list(obj = obj,
                      as_sparse = as_sparse,
                      verbose = verbose)
  #### Convert matrix ####
  if(isTRUE(as_delayedarray)) {
    obj1$data <- to_delayedarray(obj = obj1$data,
                                 verbose = verbose)
  }
  #### Create se object ####
  obj2 <- SummarizedExperiment::SummarizedExperiment(
    assays = obj1$data,
    colData = obj1$obs,
    rowData = obj1$var
  )
  obj2 <- check_se_rownames(obj2, verbose = verbose)
  return(obj2)
}
