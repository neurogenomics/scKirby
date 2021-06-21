
is.EWCElist <- function(object){
  class(object)[1]=="list" & all(c("exp","annot") %in% names(object))
}
