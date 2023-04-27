#' Set unstructured data
#'
#' Set the sample observations slot (i.e. cell metadata)
#'  in any single-cell object that has one.
#' @param obs Observation metadata
#' as a \link[base]{data.frame} with samples as row names.
#' @returns Single-cell object.
#'
#' @export
#' @examples
#' obj <- example_obj("ad")
#' obs <- get_obs(obj)
#' obs$new_col <- c(1,2)
#' obj2 <- set_obs(obj = obj,
#'                 obs = obs)
set_obs <- function(obj,
                    obs,
                    verbose=TRUE){

  messager("Setting cell metadata (obs) in obj.",v=verbose)
  obs <- check_metadata_rownames(d = obs,
                                 verbose = verbose)
  #### anndata ####
  if(is_class(obj,"anndata")){
    obj$obs <- obs[obj$obs_names,]
  #### seurat ####
  } else if(is_class(obj,"seurat")){
    obj@meta.data <- obs[colnames(obj),]
  }
  return(obj)
}
