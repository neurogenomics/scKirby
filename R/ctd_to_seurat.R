#' Convert: \code{CellTypeDataSet} ==> \code{Seurat}
#'
#' @inheritParams converters
#' @export
#' @examples
#' obj <- example_obj("ctd")
#' obj2 <- ctd_to_seurat(obj)
ctd_to_seurat <- function(obj,
                          verbose=TRUE){
  obj <- ctd_to_se(obj = obj,
                   verbose = verbose)
  obj <- se_to_seurat(obj = obj,
                      verbose = verbose)
  return(obj)
}
