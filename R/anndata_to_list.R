#' Convert: \code{AnnData} ==> \code{list}
#'
#' @export
#' @examples
#' obj <- example_obj("anndata")
#' obj2 <- anndata_to_list(obj)
anndata_to_list <- function(obj,
                            verbose=TRUE){
  messager("+ AnnData ==> list",v=verbose)
  list(data = Matrix::t(obj$X),
       obs = obj$obs,
       var = obj$var,
       var_features = obj$varm,
       reductions = obj$obsm,
       uns = obj$uns)
}
