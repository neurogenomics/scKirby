simplify_list <- function(l){
  if(is.list(l) && length(l)==1) l <- l[[1]]
  return(l)
}
