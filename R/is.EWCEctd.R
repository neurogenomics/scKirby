

is.EWCEctd <- function(object){
  tryCatch(expr = {
    (class(object)[1]=="list" & length(object)>1) &
      all(c("mean_exp","specificity") %in% names(object[[1]]))
  }, error=function(x)return(F))
}
