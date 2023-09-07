#' Set graphs
#'
#' Set graph objects in any single-cell object class.
#' @param g A graph.
#' @param key Name of the graph key to set \code{g} as within \code{obj}.
#' @inheritParams converters
#' @returns Named list of graph objects.
#'
#' @export
#' @importFrom methods is
#' @examples
#' obj <- example_obj("seurat")
#' g <- get_cor(obj)
#' obj2 <- set_graph(obj=obj, g=g, key="cor_graph")
set_graph <- function(obj,
                      g,
                      key,
                      verbose = TRUE) {
  force(obj);force(g);force(key);
  #### Check graph ####
  g <- to_graph(obj = g,
                verbose = verbose)
  #### Set graph ####
  if (is_class(obj,"seurat")) {
    messager("Setting graph in Seurat obj:",key,v=verbose)
    obj@graphs[[key]] <- g
  } else {
    stopper("obj type not suppported yet.")
  }
  return(obj)
}
