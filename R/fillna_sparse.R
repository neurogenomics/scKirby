fillna_sparse <- function(X,
                          drop_allna=TRUE){
  if(any(is.na(X))){
    ## Drop cols with only NAs
    if(isTRUE(drop_allna)){
      all_na <- apply(X,2,function(x){all(is.na(x))})
      if(sum(all_na)>0){
        X <- X[,!all_na,drop=FALSE]
      }
    }
    ## Fill cols with some NAs
    X[is.na(X)] <- 0
  }
  return(X)
}
