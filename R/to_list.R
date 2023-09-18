#' Convert: to \code{list}
#'
#' @describeIn converters
#' Convert any object to \link[base]{list} format.
#' @inheritParams converters
#' @returns A \link[base]{list} object.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- to_list(obj)
to_list <- function(obj,
                    verbose=TRUE){
  #### Template ####
  # list(data = NULL,
  #      obs = NULL,
  #      var = NULL,
  #      obsm = NULL,
  #      varm = NULL,
  #      graphs = NULL,
  #      uns = NULL)

  #### matrix/data.frame ####
  if(is_class(obj,"matrix")){
    l <- list(data=to_sparse(obj = obj,
                             verbose = verbose))
  #### loom ####
  } else if(is_class(obj,"loom")){
    l <- loom_to_list(obj = obj,
                      verbose = verbose)
  #### SummarizedExperiment ####
  } else if(is_class(obj,"se")){
    l <- se_to_list(obj = obj,
                    verbose = verbose)
  #### Seurat ####
  } else if(is_class(obj,"seurat")){
    l <- seurat_to_list(obj = obj,
                        verbose = verbose)
  #### anndata ####
  } else if(is_class(obj,"anndata")){
    l <- anndata_to_list(obj = obj,
                         verbose = verbose)
  #### CellDataSet ####
  } else if(is_class(obj,"cds")){
    l <- cds_to_list(obj = obj,
                     verbose = verbose)
  #### OTHER ####
  } else {
    obj <- to_se(obj = obj,
                 verbose = verbose)
    l <- se_to_list(obj = obj,
                    verbose = verbose)
  }
  return(l)
}
