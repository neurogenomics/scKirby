#' Convert: \code{anndata} ==> \code{list}
#'
#' @export
#' @examples
#' obj <- example_obj("anndata")
#' obj2 <- anndata_to_list(obj)
anndata_to_list <- function(obj){
  list(data = Matrix::t(obj$X),
       obs = obj$obs,
       var = obj$var,
       var_features = NULL,
       reductions = NULL,
       graphs = NULL)
}
