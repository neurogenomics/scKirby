#' Get observations metadata
#'
#' Extract sample observations (i.e. cell) metadata from any single-cell object.
#' @param rownames_col Name of the column to use as row names in the metadata
#' (i.e. unique cell IDs/barcodes).
#' @inheritParams converters
#' @returns An observation (sample) metadata data.frame.
#'
#' @export
#' @examples
#' obj <- example_obj("scle")
#' obs <- get_obs(obj)
get_obs <- function(obj,
                    rownames_col=NULL,
                    verbose=TRUE){
  # devoptera::args2vars(get_obs)

  #### Matrix ####
  if(is_class(obj,"matrix")){
    obs <- data.frame(variable = colnames(obj),
                      row.names = colnames(obj))
  #### loom ####
  } else if(is_class(obj,"loom")){
    obs <- as.data.frame(obj[["col_attrs"]])
  #### SummarizedExperiment ####
  } else if(is_class(obj,"se")){
    obs <- obj@colData
  #### Seurat ####
  } else if(is_class(obj,"seurat")){
    ## Seurat V1
    if(methods::is(obj,"seurat")){
      obs <- obj@data.info
    ## Seurat V2+
    } else {
      obs <- obj@meta.data
    }
  #### h5Seurat ####
  } else if(is_class(obj,"h5seurat")){
  obs <- as.data.frame(obj[["meta.data"]])
  #### anndata ####
  } else if (methods::is(obj,"DimReduc")) {
    messager("Using embedding rownames as metadata.",v=verbose)
    obs <- data.frame(
      id = rownames(obj@cell.embeddings),
      ## Redundant but extra column prevents df
      ## from turning into list sometimes.
      label_phe_code = rownames(obj@cell.embeddings),
      row.names = rownames(obj@cell.embeddings)
    )
  } else if(is_class(obj,"anndata")){
    obs <- obj$obs
  #### CellDataSet ####
  } else if(is_class(obj,"cds")){
    obs <- obj@phenoData@data
  #### list ####
  } else if(is_class(obj,"list")){
    #### File path ####
    if(is.character(obj$obs)){
      obs <- read_data(path = obj$obs,
                       verbose = verbose,
                       as_sparse = FALSE)
    } else {
      obs <- obj$obs
    }
  #### OTHER ####
  } else {
    stopper("Unable to get `obs` from object.")
  }
  #### Add rownames ####
  obs <- check_metadata_rownames(d = obs,
                                 rownames_col = rownames_col,
                                 verbose = verbose)
  #### Return ####
  return(obs)
}
