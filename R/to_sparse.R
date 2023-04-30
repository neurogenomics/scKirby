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
                      allow_delayed_array = TRUE,
                      verbose = TRUE){

  #### Return directly ####
  is_sparse_matrix <- utils::getFromNamespace("is_sparse_matrix","orthogene")
  if(isTRUE(is_sparse_matrix(X = obj))) return(obj)
  #### Turn first column into rownames ####
  if(is.null(rownames(obj)) ||
     is.character(obj[[1]])){
    obj <- as.matrix(obj[,-1]) |> `rownames<-`(obj[[1]])
  }
  #### Convert to sparse matrix ####
  utils::getFromNamespace("to_sparse","orthogene")(
    gene_df2 = obj,
    allow_delayed_array = allow_delayed_array,
    verbose = verbose)
}
