#' Extract correlation matrix
#'
#' Extract correlation matrix stored as a graph.
#' @param assay Assay to used to compute correlation graph
#' if one does not already exist.
#' @param graph_key Name of the graph to extract.
#' @param method Pairwise correlation method.
#' @param as_graph Convert the correlation matrix to the
#' \link[SeuratObject]{Graphs} class.
#' @param return_obj Whether to return the single-cell object with a new
#' \code{graph}, or to simply return the sparse correlation matrix.
#' @inheritParams converters
#' @inheritParams calc_cor
#' @inheritParams get_obsm
#' @inheritParams SeuratObject::CreateSeuratObject
#' @inheritParams SeuratObject::Assays
#' @returns Trait-trait correlation matrix.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' Xcor <- get_cor(obj = obj,
#'                 keys = "pca")
get_cor <- function(obj,
                    keys = NULL,
                    assay = NULL,
                    slot = NULL,
                    graph_key = NULL,
                    method = "pearson",
                    return_obj = FALSE,
                    as_graph = TRUE,
                    verbose = TRUE) {
  # devoptera::args2vars(get_cor)

  if(is_class(obj,"matrix") &&
     isTRUE(return_obj)){
    messager("Cannot return object when obj is a matrix.",
             "Setting return_obj=FALSE.",v=verbose)
    return_obj <- FALSE
  }
  #### Reassign name to distinguish from other cor matrices ####
  graph_key <- infer_graph_key(obj = obj,
                               graph_key = graph_key,
                               assay = assay,
                               keys = keys,
                               ignore_has_graph = TRUE,
                               verbose = verbose)
  #### Check if graph exists ####
  Xcor <- get_graphs(obj = obj,
                    keys = graph_key,
                    verbose = verbose)
  #### if not, compute new _cor graph ####
  if(is.null(Xcor)){
    Xcor <- calc_cor(obj = obj,
                     assay = assay,
                     slot = slot,
                     keys = keys,
                     method = method)
  }
  #### Convert to sparse graph ####
  if(isTRUE(as_graph)){
    Xcor <- to_graph(Xcor)
  }
  #### Add back into object ####
  if (isTRUE(return_obj) &&
      is_class(obj,"seurat")) {
    graph_key <- infer_graph_key(obj = obj,
                                 graph_key = graph_key,
                                 assay = assay,
                                 keys = keys,
                                 ignore_has_graph = TRUE,
                                 verbose = FALSE)
    messager("Adding new graph to obj:", graph_key, v = verbose)
    obj <- set_graph(obj = obj,
                     g = Xcor,
                     key =  graph_key,
                     verbose = verbose)
    return(obj)
  } else {
    #### Return directly ####
    messager("Returning sparse correlation matrix.", v = verbose)
    return(Xcor)
  }
  return(obj2)
}
