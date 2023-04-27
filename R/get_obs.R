#' Get observations metadata
#'
#' Extract sample observations (i.e. cell) metadata from any single-cell object.
#' @param rownames_col Name of the column to use as row names in the metadata
#' (i.e. unique cell IDs/barcodes).
#' @export
#' @examples
#' obj <- example_obj("scle")
#' obs <- get_obs(obj)
get_obs <- function(obj,
                    rownames_col=NULL,
                    verbose=TRUE){
  # devoptera::args2vars(get_obs)

  #### loom ####
  if(is_class(obj,"loom")){
    obs <- as.data.frame(obj[["col_attrs"]])
  #### SummarizedExperiment ####
  } else if(is_class(obj,"se")){
    obs <- obj@colData
  #### Seurat ####
  } else if(is_class(obj,"seurat")){
    obs <- obj@meta.data
  #### h5Seurat ####
  } else if(is_class(obj,"h5seurat")){
  obs <- as.data.frame(obj[["meta.data"]])
  #### anndata ####
  } else if(is_class(obj,"anndata")){
    obs <- obj$obs
  #### CellDataSet ####
  } else if(is_class(obj,"cds")){
    obs <- obj@phenoData@data
  #### list ####
  } else if(is_class(obj,"list")){
    obs <- obj$obs
  #### OTHER ####
  } else {
    stopper("Unable to get metadata from object.")
  }
  #### Add rownames ####
  obs <- check_metadata_rownames(d = obs,
                                 rownames_col = rownames_col,
                                 verbose = verbose)
  #### Return ####
  return(obs)
}
