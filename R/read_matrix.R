read_matrix <- function(path,
                        transpose=FALSE,
                        as_sparse=TRUE,
                        verbose,
                        ...){
  if(is.character(path)){
    messager("+ Expression matrix (.mtx/.csv/.tsv) detected.",
             "Importing as sparse matrix.",v=verbose)
    obj <- data.table::fread(path,
                             ...)
  } else {
    obj <- path
  }
  #### Convert to sparse matrix ####
  if(isTRUE(as_sparse)){
    obj <- to_sparse(obj = obj,
                     verbose = verbose)
  }
  if(isTRUE(transpose)){
    messager("Transposing matrix.",v=verbose)
    obj <- Matrix::t(obj)
  }
  return(obj)
}
