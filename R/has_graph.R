#' Has graph
#'
#' Check whether an object has a particular graph.
#' @param graph_keys The names of specific graphs to extract.
#' @inheritParams converters
has_graph <- function(obj,
                      graph_keys,
                      verbose = TRUE){

    g <- get_obsp(obj = obj,
                    verbose = verbose)
    graph_keys %in% names(g)
}
