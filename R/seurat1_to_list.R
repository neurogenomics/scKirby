#' Convert: \code{Seurat (V1)} ==> \code{list}
#'
#' @export
seurat1_to_list <- function(obj,
                            verbose=TRUE){

  messager("+ Seurat (V1) ==> list",v=verbose)
  list(data = get_data(obj = obj, verbose = verbose),
       obs = get_obs(obj = obj, verbose = verbose),
       var = get_var(obj = obj, verbose = verbose),
       reductions = get_reductions(obj = obj, verbose = verbose),
       graphs = get_graphs(obj = obj, verbose = verbose)
  )
}
