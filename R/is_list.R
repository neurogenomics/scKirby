is_list <- function(obj,
                    validate=FALSE){
  is.list(obj) &&
    !is.data.frame(obj) &&
    !is_sparse_matrix(obj) &&
    !methods::is(obj,"matrix") &&
    if(isTRUE(validate)){
      all(c("data","obs") %in% names(obj))
    } else {
      TRUE
    }
}
