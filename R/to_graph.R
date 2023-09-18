#' To graph
#'
#' @describeIn converters
#' Convert a matrix to a Graph object (using \link[SeuratObject]{as.Graph}).
#' @inheritParams converters
#' @inheritDotParams get_x
#' @returns A sparse \link[SeuratObject]{as.Graph} object.
#'
#' @export
#' @importFrom SeuratObject as.Graph
#' @examples
#' obj <- example_obj("matrix")
#' obj2 <- to_graph(obj)
to_graph <- function(obj,
                     verbose = TRUE,
                     ...){
  Xl <- get_x(obj = obj,
              verbose = verbose,
              ...)
  lapply(Xl,Seurat::as.Graph)
}
