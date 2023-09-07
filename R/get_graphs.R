#' Get graphs
#'
#' Get graph objects from any single-cell object class.
#' @param keys The keys of reductions to extract from.
#' @inheritParams converters
#' @inheritParams get_n_elements
#' @returns Named list of graph objects.
#'
#' @export
#' @importFrom methods is
#' @examples
#' obj <- example_obj("seurat")
#' g <- get_graphs(obj)
get_graphs <- function(obj,
                       keys = NULL,
                       n = NULL,
                       verbose = TRUE) {

  if (is_class(obj,"seurat")) {
    ## Seurat V1
    if(methods::is(obj,"seurat")){
      g <- list(snn.sparse=obj@snn.sparse)
    ## Seurat V2+
    } else {
      all_keys <- rev(names(obj@graphs))
      if(!is.null(keys)){
        all_keys <- all_keys[tolower(all_keys) %in% tolower(keys)]
      }
      if(length(all_keys)==0){
        messager("None of the requested graphs can be found.",
                 "Returning NULL.",v=verbose)
        return(NULL)
      }
      g <- obj@graphs[all_keys]
    }

  } else if (methods::is(obj, "Graph")) {
    messager("Using obj as graph.", v = verbose)
    g <- obj
  } else {
    messager("No graph found. Returning NULL.", v = verbose)
    return(NULL)
  }
  #### Return as a named list (1 per assay), unless there's only 1 assay ####
  g <- get_n_elements(l = g,
                      n = n,
                      verbose = verbose)
  return(g)
}
