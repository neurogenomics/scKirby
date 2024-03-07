#' @describeIn anndata_to anndata_to
#'
#' @export
#' @examples
#' obj <- example_obj("anndata")
#' obj2 <- anndata_to_list(obj)
anndata_to_list <- function(obj,
                            verbose=TRUE){

  messager("+ AnnData ==> list",v=verbose)
  list(data = get_x(obj = obj,
                       verbose = verbose),
       obs = get_obs(obj = obj,
                     verbose = verbose),
       var = get_var(obj = obj,
                     verbose = verbose),
       obsm = get_obsm(obj = obj,
                             verbose = verbose),
       varm = get_varm(obj = obj,
                       verbose = verbose),
       uns = obj$uns)
}
