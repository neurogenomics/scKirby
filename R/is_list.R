is_list <- function(obj){
  is.list(obj) &&
    all(c("data","obs") %in% names(obj))
}
