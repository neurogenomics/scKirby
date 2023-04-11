#' Convert: \code{CellDataSet} ==> \code{list}
#'
#' @export
#' @examples
#' obj <- example_obj("cds")
#' obj2 <- cds_to_list(obj)
cds_to_list <- function(obj){
  list(data=Biobase::exprs(obj),
       obs=Biobase::pData(obj),
       var=Biobase::fData(obj))
}
