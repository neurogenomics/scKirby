read_matrix <- function(path,
                        verbose,
                        ...){
  messager("+ Expression matrix (.mtx) detected.",
           "Importing as sparse matrix.",v=verbose)
  obj <- data.table::fread(path,
                           verbose = verbose,
                           ...)
  as.matrix(obj[,-1], obj[[1]]) |>
    Matrix::Matrix(sparse=TRUE)
}
