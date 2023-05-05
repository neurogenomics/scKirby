#' To sparse matrix
#'
#' Convert a matrix or data.frame to a sparse matrix.
#' @returns Sparse matrix
#'
#' @export
#' @examples
#' obj <- data.frame(gene=paste0("gene",seq_len(10)),
#'                   matrix(data = 0,nrow=10,ncol = 20))
#' obj2 <- to_sparse(obj)
to_sparse <- function(obj,
                      transpose = FALSE,
                      allow_delayed_array = TRUE,
                      verbose = TRUE){

  #### Check if already a sparse matrix ####
  is_sparse_matrix <- utils::getFromNamespace("is_sparse_matrix","orthogene")
  if(isFALSE(is_sparse_matrix(obj))){
    #### Turn first column into rownames ####
    if(is.null(rownames(obj)) ||
       is.character(obj[[1]])){
      obj <- as.matrix(obj[,-1]) |> `rownames<-`(obj[[1]])
    }
    #### Convert to sparse matrix ####
    obj <- utils::getFromNamespace("to_sparse","orthogene")(
      gene_df2 = obj,
      allow_delayed_array = allow_delayed_array,
      verbose = verbose)
  }
  ##### Transpose ####
  if(isTRUE(transpose)){
    messager("Transposing matrix.",v=verbose)
    obj <- Matrix::t(obj)
  }
  #### Return ####
  return(obj)
}
