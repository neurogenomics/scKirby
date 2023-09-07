#' Convert: \code{SummarizedExperiment} ==> \code{list}
#'
#' @inheritParams converters
#' @export
#' @examples
#' obj <- example_obj("se")
#' obj2 <- se_to_list(obj)
se_to_list <- function(obj,
                       verbose=TRUE){

  messager("+ AnnData ==> Seurat",v=verbose)
  list(data = get_x(obj, verbose = verbose),
       obs = get_obs(obj, verbose = verbose),
       var = get_var(obj, verbose = verbose))
}
