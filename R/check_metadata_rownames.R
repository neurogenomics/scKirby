check_metadata_rownames <- function(d,
                                    rownames_col=NULL,
                                    verbose=TRUE){
  #### Recursion ####
  if(is.list(d) &&
     !is_class(d,"matrix")){
    d <- lapply(d, check_metadata_rownames, verbose=verbose)
    return(d)
  }
  if(is.null(rownames(d))){
    valid_rn <- sapply(d,function(x){length(unique(x))==nrow(d)})
    if(sum(valid_rn)>0){
      if(rownames_col %in% valid_rn){
        rn_col <- rownames_col
      } else {
        rn_col <- names(valid_rn)[[1]]
      }
      messager("Assigning column",shQuote(rn_col),
               "as metadata row names.",v=verbose)
      rownames(d) <- d[[rn_col]]
    } else {
      wrn <- paste(
        "No columns could be identified as metadata row names.",
        "Please assign this manually."
      )
      warning(wrn)
    }
  }
  return(d)
}
