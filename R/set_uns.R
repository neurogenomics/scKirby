#' Set unstructured data
#'
#' Set the unstructured data layer in any single-cell object that has one.
#' @param uns Unstructured data to be stored in object.
#' @param key Name of the list element to store \code{uns} in.
#' @inheritParams converters
#' @returns Single-cell object.
#'
#' @export
#' @examples
#' obj <- example_obj("ad")
#' uns <- list("extra_info"=mtcars)
#' obj2 <- set_uns(obj = obj,
#'                 uns = mtcars,
#'                 key = "extra_info")
set_uns <- function(obj,
                    uns,
                    key,
                    verbose=TRUE){

  messager("Setting unstructured data in obj.",v=verbose)
  if(is_class(obj,"anndata")){
    obj$uns[[key]] <- uns
  } else if(is_class(obj,"seurat")){
    obj@misc[[key]] <- uns
  }
  return(obj)
}
