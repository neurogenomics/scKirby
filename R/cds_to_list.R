#' Convert: \code{CellDataSet} ==> \code{list}
#'
#' @export
#' @examples
#' obj <- example_obj("cds")
#' obj2 <- cds_to_list(obj)
cds_to_list <- function(obj,
                        verbose = TRUE){
  messager("+ CellDataSet ==> list",v=verbose)
  # Biobase::ExpressionSet()
  list(data = get_data(obj, verbose = verbose),
       obs = get_obs(obj, verbose = verbose),
       var = get_var(obj, verbose = verbose),
       obsm = get_obsm(obj,
                       verbose = verbose)
       )
}
