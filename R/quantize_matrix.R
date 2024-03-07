#' Quantize matrix
#'
#' @description Quantize matrix
#' @inheritParams quantize_vector
#' @inheritParams to_sparse
#'
#' @export
#' @examples
#' X <- example_obj("matrix")
#' Xq <- quantize_matrix(X)
quantize_matrix <- function(X,
                            n=40,
                            as_sparse=TRUE,
                            verbose=TRUE){
    ## Quantize matrix
    Xq <- apply(X, 2,
                FUN = quantize_vector,
                n = n) |>
        as.matrix() |>
        `row.names<-`(rownames(X))
    ## Convert to sparse matrix
    if(isTRUE(as_sparse)){
        Xq <- to_sparse(
            obj = Xq,
            verbose = verbose
        )
    }
    ## Return
    return(Xq)
}
