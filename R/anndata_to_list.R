#' Convert: \code{AnnData} ==> \code{list}
#'
#' @export
#' @examples
#' obj <- example_obj("anndata")
#' obj2 <- anndata_to_list(obj)
anndata_to_list <- function(obj,
                            verbose=TRUE){
  messager("+ AnnData ==> list",v=verbose)
  list(data = get_data(obj = obj,
                       verbose = verbose),
       obs = get_obs(obj = obj,
                     verbose = verbose),
       var = get_var(obj = obj,
                     verbose = verbose),
       var_features = obj$varm,
       reductions = get_reductions(obj = obj,
                                   verbose = verbose),
       uns = obj$uns)
}
