#' To \code{AnnData}
#'
#' @describeIn converters
#' Convert any single-cell object to \link[anndata]{AnnData} format.
#' @param reimport Save and re-import the \link[anndata]{AnnData} object
#' into R to ensure all data has been converted from Python-native to
#' R-native objects
#' (e.g. pandas data.frames vs. R \link[base]{data.frames}).
#' @inheritParams converters
#' @returns An \link[anndata]{AnnData} object.
#'
#' @export
#' @examples
#' obj <- example_obj("cds")
#' obj2 <- to_anndata(obj)
to_anndata <- function(obj,
                       save_path = tempfile(fileext = ".h5ad"),
                       reimport = FALSE,
                       verbose = TRUE){

  #### Check if class is supported ####
  check_supported(obj)
  #### seurat ####
  if(is_class(obj,"seurat")){
    obj2 <- seurat_to_anndata(obj = obj,
                              save_path = save_path,
                              reimport = reimport,
                              verbose = verbose)
  #### h5seurat  ####
  } else if(is_class(obj,"h5seurat")){
    obj2 <- h5seurat_to_anndata(obj = obj,
                                save_path = save_path,
                                reimport = reimport,
                                verbose = verbose)
    #### CTD  ####
  } else if(is_class(obj,"ctd")){
    # obj2 <- ctd_to_anndata(obj = obj,
    #                       verbose = verbose)
  #### Matrices ####
  } else if(is_class(obj,"matrix")){
    # obj2 <- matrix_to_anndata(obj = obj,
    #                           verbose = verbose)
  #### list ####
  } else if(is_class(obj,"list")){
    obj2 <- list_to_anndata(obj = obj,
                            verbose = verbose)
  #### anndata ####
  } else if(is_class(obj,"anndata")){
    messager("obj already in AnnData format.")
    obj2 <- obj
  #### loom ####
  } else if(is_class(obj,"loom")){
    # obj2 <- loom_to_anndata(obj = obj,
    #                         verbose = verbose)
  #### SummarizedExperiment ####
  } else if(is_class(obj,"se")){
    obj2 <- se_to_anndata(obj = obj,
                          save_path = save_path,
                          reimport = reimport,
                          verbose = verbose)
  #### CellDataSet ####
  } else if(is_class(obj,"cds")){
    obj2 <- cds_to_anndata(obj = obj,
                           save_path = save_path,
                           reimport = reimport,
                           verbose = verbose)
  #### OTHER ####
  } else {
    l <- to_list(obj = obj,
                 verbose = verbose)
    obj2 <- list_to_anndata(obj = l,
                            verbose = verbose)
  }
  #### Return ####
  return(obj2)
}
