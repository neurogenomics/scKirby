is_list <- function(obj){
  is.list(obj) &&
    !is.data.frame(obj) &&
    all(c("data","obs") %in% names(obj))
}
