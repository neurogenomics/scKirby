#' To graph
#'
#' Convert a matrix to a Graph object (using \link[Seurat]{as.Graph}).
#' @inheritParams converters
#' @returns A sparse \link[Seurat]{as.Graph} object.
#'
#' @export
#' @importFrom Seurat as.Graph
#' @examples
#' obj <- example_obj("matrix")
#' obj2 <- to_graph(obj)
to_graph <- function(obj,
                     verbose = TRUE){
  Seurat::as.Graph(obj)
}
