#' To graph
#'
#' @describeIn converters
#' Convert a matrix to a Graph object (using \link[SeuratObject]{as.Graph}).
#' @param Xl A named list of matrices.
#' @inheritParams converters
#' @inheritDotParams get_x
#' @returns A sparse \link[SeuratObject]{as.Graph} object.
#'
#' @export
#' @examples
#' obj <- example_obj("matrix")
#' obj2 <- to_graph(obj)
to_graph <- function(Xl,
                     verbose = TRUE,
                     ...){
  messager_to()
  Xl <- to_list(Xl)
  lapply(Xl,Seurat::as.Graph)
}
