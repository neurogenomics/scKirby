check_se_rownames <- function(obj,
                               rownames_var=NULL,
                               remove_duplicates=TRUE,
                               verbose=TRUE){

  messager("+ Checking obj rownames.",v=verbose)
  rowDat <- SummarizedExperiment::rowData(obj)
  rownames_var <- if(is.null(rownames_var)) "Gene" else rownames_var
  if(is.null( S4Vectors::rownames(rowDat)) ){
    if(rownames_var %in% colnames(rowDat)){
      messager("+ Assigning rownames to:",rownames_var,v=verbose)
      # IMPORTANT! only S4Vectors::rownames can assign row names.
      # Automatically assigns same rownames to assay and rowData
      S4Vectors::rownames(obj) <- rowDat[[rownames_var]]
    }else {
      wrn <- paste(
        "Cannot identify rownames.",
        "Please set rownames (rowDat) first with:",
        "S4Vectors::rownames(obj) <- gene_names")
      warning(wrn)
    }
  }
  if(isTRUE(remove_duplicates)){
    dup_rows <- base::duplicated(S4Vectors::rownames(obj))
    if(sum(dup_rows)>0){
      messager("+ Removing duplicate gene rows.",v=verbose)
      obj <- obj[!dup_rows,]
    }
  }
  return(obj)
}
