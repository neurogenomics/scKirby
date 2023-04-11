#' Convert: \code{CellTypeDataSet} ==> \code{Seurat}
#'
#' @export
#' @examples
#' data("example_ctd")
#' obj2 <- ctd_to_seurat(example_ctd)
ctd_to_seurat <- function(obj,
                          verbose=TRUE){
  obj <- ctd_to_sce(obj = obj,
                    verbose = verbose)
  obj <- sce_to_seurat(obj = obj,
                       verbose = verbose)
  return(obj)
}
