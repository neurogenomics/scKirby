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
#' g <- get_obsp(obj)
get_obsp <- function(obj,
                     keys = NULL,
                     n = NULL,
                     as_graph = FALSE,
                     verbose = TRUE) {
  if(is_class(obj,"matrix")){
    obsp <- list(obsp=obj)
  } else if(is_class(obj,"matrix_list")){
    obsp <- obj
  } else if (methods::is(obj, "Graph")) {
    messager("Using obj as graph.", v = verbose)
    obsp <- list("graph"=obj)
  } else if (is_class(obj,"list")){
    obsp <- obj$obsp
  } else if (is_class(obj,"seurat")) {
    ## Seurat V1
    if(methods::is(obj,"seurat")){
      obsp <- list(snn.sparse=obj@snn.sparse)
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
      obsp <- obj@graphs[all_keys]
    }
  } else if(is_class(obj,"anndata")){
    obsp <- obj$obsp
  }else {
    messager("No graph found. Returning NULL.", v = verbose)
    return(NULL)
  }
  ### Convert to graphs
  if(isTRUE(as_graph)){
    obsp <- to_graph(obsp)
  }
  #### Return as a named list (1 per assay), unless there's only 1 assay ####
  obsp <- get_n_elements(l = obsp,
                         n = n,
                         verbose = verbose)
  return(obsp)
}
