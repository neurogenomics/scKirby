#' Set unstructured data
#'
#' Set the feature variable layer (i.e. gene metadata)
#' in any single-cell object that has one.
#' @param var Variable metadata
#' as a \link[base]{data.frame} with features as row names.
#' @inheritParams converters
#' @returns Single-cell object.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' var <- get_var(obj)
#' var$new_col <- c(1,2)
#' obj2 <- set_var(obj = obj,
#'                 var = var)
set_var <- function(obj,
                    var,
                    verbose=TRUE){

  messager("Setting feature variable metadata (var) in obj.",v=verbose)
  var <- check_metadata_rownames(d = var,
                                 verbose = verbose)
  #### anndata ####
  if(is_class(var,"anndata")){
    obj$var <- var[obj$var_names,]
  #### seurat ####
  } else if(is_class(obj,"seurat")){
    for(a in names(obj@assays)){
      if(a %in% names(var)){
        obj@assays[[a]]@meta.features <- var[[a]][rownames(obj),,drop=FALSE]
      }
    }
  }
  return(obj)
}
