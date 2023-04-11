check_se_colnames <- function(obj,
                              colnames_var=NULL,
                              verbose=TRUE){
  messager("Checking obj colnames.",v=verbose)
  colDat <- SummarizedExperiment::colData(obj)
  if(!is.null(colnames_var)){
    if(colnames_var %in% colnames(colDat)){
      messager("Assigning colnames to:", colnames_var,v=verbose)
      S4Vectors::colnames(obj) <- colDat[[colnames_var]]
    } else {
      messager(colnames_var,"not found in colData.",v=verbose)
    }
  }
  return(obj)
}



