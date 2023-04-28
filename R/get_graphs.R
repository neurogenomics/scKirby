#' Get graphs
#'
#' Get graph objects from any single-cell object class.
#' @param names The names of reductions to extract from.
#' @param verbose Print messages.
#'
#' @export
#' @importFrom methods is
#' @examples
#' obj <- example_obj("seurat")
#' g <- get_graphs(obj)
get_graphs <- function(obj,
                       names = NULL,
                       verbose = TRUE) {

  if (is_class(obj,"seurat")) {
    ## Seurat V1
    if(methods::is(obj,"seurat")){
      g <- list(snn.sparse=obj@snn.sparse)
    ## Seurat V2+
    } else {
      all_names <- rev(names(obj@graphs))
      if(!is.null(names)){
        all_names <- all_names[tolower(all_names) %in% tolower(names)]
      }
      g <- obj@graphs[all_names]
    }

  } else if (methods::is(obj, "Graph")) {
    messager("Using obj as graph.", v = verbose)
    g <- obj
  } else {
    messager("No graph found. Returning NULL.", v = verbose)
    g <- NULL
  }
  return(g)
}
