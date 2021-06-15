check_sce_rownames <- function(sce,
                               rownames_var=NULL,
                               remove_duplicates=T,
                               verbose=T){
  printer("Checking SCE rownames...",v=verbose)
  rowDat <- SummarizedExperiment::rowData(sce)
  rownames_var <- if(is.null(rownames_var)) "Gene" else rownames_var
  if( is.null( S4Vectors::rownames(rowDat)) ){
    if(rownames_var %in% colnames(rowDat)){
      printer("+ Assigning rownames to:",rownames_var,v=verbose)
      # IMPORTANT! only S4Vectors::rownames can assign row names.
      S4Vectors::rownames(sce) <- rowDat[[rownames_var]] # Automatically assigns same rownames to assay and rowData
    }else {
      warning("Cannot identify rownames. Please set rownames (rowDat) first with: S4Vectors::rownames(sce) <- gene_names")
    }
  }
  if(remove_duplicates){
    printer("+ Removing duplicate gene rows.",v=verbose)
    sce <- sce[!base::duplicated(S4Vectors::rownames(sce)),]
  }
  return(sce)
}
