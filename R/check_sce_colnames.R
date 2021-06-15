
check_sce_colnames <- function(sce,
                               colnames_var=NULL,
                               verbose=T){
  printer("Checking SCE colnames.",v=verbose)
  colDat <- SummarizedExperiment::colData(sce)
  if(!is.null(colnames_var)){
    if(colnames_var %in% colnames(colDat)){
      printer("Assigning colnames to:", colnames_var,v=verbose)
      S4Vectors::colnames(sce) <- colDat[[colnames_var]]
    } else {
      printer(colnames_var,"not found in colData.",v=verbose)
    }
  }
  return(sce)
}



