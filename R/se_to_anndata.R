#' Convert: \code{SummarizedExperiment} ==> \code{AnnData}
#'
#' @export
#' @import sceasy
#' @examples
#' obj <- example_obj("se")
#' obj2 <- se_to_anndata(obj)
se_to_anndata <- function(obj,
                          verbose=TRUE){
  messager("+ SummarizedExperiment ==> AnnData",v=verbose)
  obj <- se_to_sce(obj = obj,
                   verbose = verbose)
  sceasy::convertFormat(obj, from = "sce", to="anndata")
}


