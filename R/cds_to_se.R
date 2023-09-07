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

  messager("+ CellDataSet ==> SummarizedExperiment",v=verbose)
  X <- Biobase::exprs(obj)
  obs <- Biobase::pData(obj)
  var <- Biobase::fData(obj)
  #### Convert matrix ####
  if(isTRUE(as_sparse)) X <- to_sparse(obj = X, verbose = verbose)
  if(isTRUE(as_delayedarray)) X <- to_delayedarray(obj = X , verbose = verbose)
  #### Create se object ####
  obj2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(raw = X),
    colData = obs,
    rowData = var
  )
  obj2 <- check_se_rownames(obj2, verbose = verbose)
  return(obj2)
}
