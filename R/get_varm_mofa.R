get_varm_mofa <- function(obj,
                          verbose=TRUE){

  varm <- obj@expectations$W
  if(is.list(varm) && length(varm)>1){
    messager(">1 group found in MOFA loadings.",
             "Using only the first:",shQuote(names(varm)[1]),v=verbose)
    varm <- varm[[1]]
  }
  varm <- list(mofa=as.matrix(varm))
  return(varm)
}
