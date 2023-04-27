#' Convert: \code{Seurat} ==> \code{list}
#'
#' @export
#' @importFrom methods .hasSlot slot
#' @importFrom stats setNames
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- seurat_to_list(obj)
seurat_to_list <- function(obj,
                           verbose=TRUE){

 messager("+ Seurat ==> list",v=verbose)
 list(
   data = get_data(obj = obj,
                   verbose = verbose),
   obs = get_obs(obj = obj,
                 verbose = verbose),
   var = get_var(obj = obj,
                 verbose = verbose),
   var_features = lapply(obj@assays,function(a){a@var.features}),
   reductions = get_reductions(obj = obj,
                               verbose = verbose),
   graphs = get_graphs(obj = obj,
                       verbose = verbose)
 )
}
