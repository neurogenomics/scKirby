#' Convert: to \code{list}
#'
#' Convert any object to \code{list} format.
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' l <- to_list(obj)
to_list <- function(obj,
                    verbose=TRUE){
  #### Template ####
  # list(assays = NULL,
  #      obs = NULL,
  #      var = NULL,
  #      var_features = NULL,
  #      reductions = NULL,
  #      graphs = NULL)

  #### matrix/data.frame ####
  if(is_class(obj,"matrix")){
    l <- list(data=obj)
  #### loom ####
  } else if(is_class(obj,"loom")){
    l <- loom_to_list(obj = obj,
                      verbose = verbose)
  #### SummarizedExperiment ####
  } else if(is_class(obj,"se")){
    l <- se_to_list(obj = obj)
  #### Seurat ####
  } else if(is_class(obj,"seurat")){
    l <- seurat_to_list(obj = obj)
  #### anndata ####
  } else if(is_class(obj,"anndata")){
    l <- anndata_to_list(obj = obj)
  #### OTHER ####
  } else if(is_class(obj,"cds")){
    l <- cds_to_list(obj = obj)
  } else {
    obj <- ingest_data(obj = obj,
                       verbose = verbose)
    l <- se_to_list(obj = obj)
  }
  return(l)
}
