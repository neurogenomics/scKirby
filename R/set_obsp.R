#' Set graphs
#'
#' Set graph objects in any single-cell object class.
#' @param obsp A observation x observation graph.
#' @param key Name of the graph key to set \code{obsp} as within \code{obj}.
#' @inheritParams converters
#' @returns Named list of graph objects.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' g <- get_cor(obj)
#' obj2 <- set_obsp(obj=obj, obsp=obsp, key="cor_graph")
set_obsp <- function(obj,
                     obsp,
                     key,
                     verbose = TRUE) {
  force(obj);force(obsp);force(key);
  #### Check graph ####

  #### Set graph ####
  if (is_class(obj,"seurat")) {
    messager("Setting graph in Seurat obj:",key,v=verbose)
    obsp <- to_graph(obj = obsp,
                     verbose = verbose)
    obj@graphs[[key]] <- obsp
  } else if (is_class(obj,"anndata")) {
    messager("Setting graph in AnnData obj:",key,v=verbose)
    obj$obsp[[key]] <- obsp
  } else {
    stopper("obj type not suppported yet.")
  }
  return(obj)
}
