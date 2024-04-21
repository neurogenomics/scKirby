#' Is class
#'
#' When \code{group} is not supplied (default),
#' determines which group of single-cell classes the object belongs to (if any).
#' When \code{group} is supplied,
#' determine whether a single-cell object belongs to a particular group of
#' object classes.
#' See \link{dict_class} for details.
#' @inheritParams converters
#' @inheritParams dict_class
#'
#' @export
#' @examples
#' obj <- example_obj("Seurat")
#' is_class(obj,"seurat")
#'
#' X <- example_obj("matrix")
#' obj <- lapply(seq(3), function(...){X})
#' is_class(obj, "matrix_list")
is_class <- function(obj,
                     group = NULL){
  # devoptera::args2vars(is_class)

  if(is.null(obj) || is.character(obj)) return(FALSE)
  cdict <- dict_class()
  matches <- lapply(cdict, function(y){
    any(sapply(y, function(z){methods::is(obj,z)}))
  })
  groups <- names(matches[unlist(matches)])
  if(is_ctd(obj)){
    groups <- c(groups,"ewce")
  }
  if(is_list(obj, validate = TRUE)){
    groups <- c(groups,"list")
  }
  if(is.list(obj) &&
     all(mapply(obj,FUN=is_class,"matrix"))){
    groups <- c(groups,"matrix_list")
  }
  if(is.null(group)){
    return(groups)
  } else {
    #### Exceptions ####
    # if(all(c("loom","h5seurat") %in% groups)) groups <- "h5seurat"
    return(tolower(group) %in% groups)
  }
}
