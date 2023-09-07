#' Convert: \code{Seurat (V1)} ==> \code{list}
#'
#' @inheritParams converters
#' @export
seurat1_to_list <- function(obj,
                            verbose=TRUE){

  messager("+ Seurat (V1) ==> list",v=verbose)
  list(data = get_x(obj = obj, verbose = verbose),
       obs = get_obs(obj = obj, verbose = verbose),
       var = get_var(obj = obj, verbose = verbose),
       obsm = get_obsm(obj = obj, verbose = verbose),
       varm = get_varm(obj = obj, verbose = verbose),
       graphs = get_graphs(obj = obj, verbose = verbose)
  )
}
