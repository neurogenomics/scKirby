#' Get assay names
#'
#' Extract assay names from any single-cell object.
#' @returns A character vector of assay names.
#' @inheritParams converters
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' nms <- get_assay_names(obj)
get_assay_names <- function(obj,
                            verbose=TRUE,
                            ...){

  if(is_class(obj,"list") ||
     is_class(obj,"matrix_list")){
    X <- get_x(obj = obj,
               verbose = verbose,
               ...)
    nms <- unique(stringr::str_split(names(X),"\\.",
                                     n = 2,
                                     simplify = TRUE)[,1])
  } else if(is_class(obj,"h5seurat")){
    nms <- names(obj[["assays"]])
  } else if(is_class(obj,"seurat")){
    nms <- Seurat::Assays(obj)
  } else if(is_class(obj,"cds")){
    nms <- names(obj@assayData)
  } else if(is_class(obj,"se") ||
            is_class(obj,"sce")){
    nms <- SummarizedExperiment::assayNames(obj)
  } else if(is_class(obj,"loom")){
    nms <- names(obj[["layers"]])
  } else if(is_class(obj,"anndata")){
    nms <- unlist(names(obj$layers))
  } else {
    messager("Unable to identify assay names from obj.",
             "Returning NULL.",v=verbose)
    nms <- NULL
  }
  return(nms)
}
