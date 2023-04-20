#' Convert: \code{CellDataSet} ==> \code{list}
#'
#' @export
#' @examples
#' obj <- example_obj("cds")
#' obj2 <- cds_to_list(obj)
cds_to_list <- function(obj,
                        verbose = TRUE){
  messager("+ CellDataSet ==> list",v=verbose)
  nms <- grep("^reducedDim",methods::slotNames(obj),value = TRUE)
  # Biobase::ExpressionSet()
  reductions <- lapply(stats::setNames(nms,
                                       nms),
                       function(nm){ methods::slot(obj, nm)})
  list(data=Biobase::exprs(obj),
       obs=Biobase::pData(obj),
       var=Biobase::fData(obj),
       reductions=reductions)
}
