read_matrix <- function(path,
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
  obj <- to_sparse(obj = obj,
                   verbose = verbose)
  return(obj)
}
