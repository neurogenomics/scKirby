check_varm <- function(varm,
                       X){
  if(!is.null(varm)){
    varm <- lapply(varm,function(x){
      #### Transpose ####
      x <- Matrix::t(x)
      #### anndata requires the dims to match ####
      if(ncol(x)!=ncol(X)){
        return(NULL)
      } else {
        return(x)
      }
    })
    varm[sapply(varm, is.null)] <- NULL
  } else {
    varm <- NULL
  }
  return(varm)
}
