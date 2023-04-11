#' Is class
#'
#' When \code{group} is not supplied (default),
#' determines which group of single-cell classes the object belongs to (if any).
#' When \code{group} is supplied,
#' determine whether a single-cell object belongs to a particular group of
#' object classes.
#' See \link[scKirby]{class_dict} for details.
#' @param obj Data object.
#' @inheritParams class_dict
#'
#' @export
#' @examples
#' obj <- example_obj("Seurat")
#' is_class(obj,"seurat")
is_class <- function(obj,
                     group = NULL){
  if(is.null(obj) || is.character(obj)) return(FALSE)
  cdict <- class_dict()
  matches <- lapply(cdict, function(y){
    any(sapply(y, function(z){methods::is(obj,z)}))
  })
  groups <- names(matches[unlist(matches)])
  if(is.null(group)){
    return(groups)
  } else {
    return(tolower(group) %in% groups)
  }
}
