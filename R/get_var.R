#' Get variable metadata
#'
#' Extract feature variable (i.e. gene) metadata from any single-cell object.
#' @param rownames_col Name of the column to use as row names in the metadata
#' (i.e. unique gene/transcript IDs).
#' @inheritParams converters
#' @inheritParams get_n_elements
#' @returns One or more variable (feature) metadata data.frames.
#'
#' @export
#' @examples
#' obj <- example_obj("cds")
#' var <- get_var(obj)
get_var <- function(obj,
                    rownames_col=NULL,
                    n=NULL,
                    simplify=TRUE,
                    verbose=TRUE){
  # devoptera::args2vars(get_var)

  #### Matrix ####
  if(is_class(obj,"matrix")){
    var <- data.frame(variable = rownames(obj),
                      row.names = rownames(obj))
  #### loom ####
  } else if(is_class(obj,"loom")){
    var <- as.data.frame(obj[["row_attrs"]])
    #### SummarizedExperiment ####
  } else if(is_class(obj,"se")){
    var <- SummarizedExperiment::rowData(obj)
    #### Seurat ####
  } else if(is_class(obj,"seurat")){
    ## Seurat V1
    if(methods::is(obj,"seurat")){
      var <- data.frame(rownames(obj@raw.data),
                        row.names = rownames(obj@raw.data))
    ## Seurat V2+
    } else {
      assays <- Seurat::Assays(obj)
      var <- lapply(stats::setNames(assays,
                                    assays),
                    function(a){
                      obj@assays[[a]]@meta.features
                    })
    }
  #### h5Seurat ####
  } else if(is_class(obj,"h5seurat")){
    assays <- obj[["assays"]]$ls()$name
    var <- lapply(stats::setNames(assays,
                                  assays),
                  function(a){
                    as.data.frame(obj[["assays"]][[a]][["meta.features"]])
                  })
    #### anndata ####
  } else if(is_class(obj,"anndata")){
    var <- obj$var
    #### CellDataSet ####
  } else if(is_class(obj,"cds")){
    var <- obj@featureData@data
    #### list ####
  } else if(is_class(obj,"list") &&
            "var" %in% names(obj)){
    #### File path ####
    if(is.character(obj$var)){
      var <- read_data(path = obj$var,
                       verbose = verbose,
                       as_sparse = FALSE)
    } else {
      var <- obj$var
    }
    #### OTHER ####
  } else {
    stopper("Unable to get `var` from object.")
  }
  #### Add rownames ####
  var <- check_metadata_rownames(d = var,
                                 rownames_col = rownames_col,
                                 verbose = verbose)
  #### Return as a named list (1 per assay), unless there's only 1 assay ####
  var <- get_n_elements(l = var,
                        n = n,
                        simplify = simplify,
                        verbose = verbose)
  #### Return ####
  return(var)
}
