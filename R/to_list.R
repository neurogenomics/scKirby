#' Convert: to \code{list}
#'
#' Convert any object to \code{list} format.
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
    obj <- ingest_data(obj = obj,
                       verbose = verbose)
    l <- se_to_list(obj = obj,
                    verbose = verbose)
  }
  return(l)
}
