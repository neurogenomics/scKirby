#' Convert: \code{loom} ==> \code{list}
#'
#' @export
#' @examples
#' library(Seurat)
#' obj <- example_obj("loom")
#' obj2 <- loom_to_list(obj)
loom_to_list <- function(obj,
                         verbose=TRUE){
  loom_to_seurat(obj,
                 verbose = verbose) |>
    seurat_to_list()
}
