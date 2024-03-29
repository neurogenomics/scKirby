#' To \code{Seurat}
#'
#' @describeIn converters
#' Convert any single-cell object to
#' \pkg{Seurat} or \link[SeuratDisk]{h5Seurat} format.
#' @param as_h5seurat Convert to the \link[SeuratDisk]{h5Seurat} class.
#' @param update Ensure the object is updated to the latest version of Seurat.
#' @inheritParams seurat_to_h5seurat
#' @inheritParams converters
#' @returns A \pkg{Seurat} or \link[SeuratDisk]{h5Seurat} object.
#'
#' @export
#' @examples
#' obj <- example_obj("list")
#' obj2 <- to_seurat(obj)
to_seurat <- function(obj,
                      as_h5seurat = FALSE,
                      update = TRUE,
                      save_path = tempfile(fileext = ".h5seurat"),
                      verbose = TRUE){
  #### Check if class is supported ####
  check_supported(obj)
  #### seurat ####
  if(is_class(obj,"seurat")){
    if(isTRUE(as_h5seurat)){
      obj2 <- seurat_to_h5seurat(obj = obj,
                                 save_path = save_path,
                                 verbose = verbose)
    } else {
      messager("+ Object already in Seurat format. Returning as-is.",
               verbose=verbose)
      return(obj)
    }
  #### h5seurat  ####
  } else if(is_class(obj,"h5seurat")){
    if(isTRUE(as_h5seurat)){
      messager("+ Object already in h5Seurat format. Returning as-is.",
               verbose=verbose)
      return(obj)
    } else {
      obj2 <- h5seurat_to_seurat(obj = obj,
                                 verbose = verbose)
    }
  #### CTD  ####
  } else if(is_class(obj,"ctd")){
    obj2 <- ctd_to_seurat(obj = obj,
                          verbose = verbose)
  #### Matrices ####
  } else if(is_class(obj,"matrix")){
    obj2 <- matrix_to_seurat(obj = obj,
                             verbose = verbose)
  #### anndata ####
  } else if(is_class(obj,"list")){
    obj2 <- list_to_seurat(obj = obj,
                           verbose = verbose)
    #### anndata ####
  } else if(is_class(obj,"anndata")){
    obj2 <- anndata_to_seurat(obj = obj,
                              verbose = verbose)
  #### loom ####
  } else if(is_class(obj,"loom")){
    obj2 <- loom_to_seurat(obj = obj,
                           verbose = verbose)
  #### SummarizedExperiment ####
  } else if(is_class(obj,"se")){
    obj2 <- se_to_seurat(obj = obj,
                         verbose = verbose)
  #### CellDataSet ####
  } else if(is_class(obj,"cds")){
    obj2 <- cds_to_seurat(obj = obj,
                          verbose = verbose)
  #### OTHER ####
  } else {
    l <- to_list(obj = obj,
                 verbose = verbose)
    obj2 <- list_to_seurat(obj = l,
                           verbose = verbose)
  }
  #### Update object ####
  if(isTRUE(update)){
    obj2 <- update_seurat(obj = obj2,
                          verbose = verbose)
  }
  #### Convert to h5seurat ####
  if(isTRUE(as_h5seurat)){
    obj2 <- seurat_to_h5seurat(obj = obj2,
                               save_path = save_path,
                               verbose = verbose)
  }
  #### Return ####
  return(obj2)
}
