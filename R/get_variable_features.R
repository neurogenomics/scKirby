#' Get variable features
#'
#' Get pre-computed variable features (e.g. genes)
#' stored within single-cell objects.
#' @param reduce Function with which to reduce variable features across assays
#'  into a single vector (e.g. \code{intersect} or \code{union}).
#' @inheritParams converters
#' @returns A named list with variable features for each assay.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' feat <- get_variable_features(obj)
get_variable_features <- function(obj,
                                  reduce=NULL){

  if(is_class(obj,"seurat")){
    feat <- lapply(obj@assays, function(x)x@var.features)
  } else {
    messager(
      "Warning:",
      "Variable features cannot be extracted from this object type yet.",
      "Returning NULL.")
    return(NULL)
  }
  #### Check if empty ####
  if(length(Reduce(union,feat))==0){
    messager("Warning:","Variable features layer was empty. Returning NULL.")
    return(NULL)
  }
  #### Reduce list ####
  if(!is.null(reduce)) {
    feat <- Reduce(reduce,feat)
    if(length(feat)==0){
      messager("Warning:","Reduced variable features had length==0.",
               "Returning NULL.")
      return(NULL)
    }
  }
  #### Return ####
  return(feat)
}
