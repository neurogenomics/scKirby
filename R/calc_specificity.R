calc_specificity <- function(X,
                             numberOfBins=40,
                             verbose=TRUE) {

  ### From EWCE::generate_celltype_data()
  messager("+ Calculating specificity matrix.",v=verbose)
  Xq <- apply(X, 2,
              FUN = EWCE::bin_columns_into_quantiles,
              numberOfBins = numberOfBins
  )
  #### Ensure rownames are kept ####
  ## Can disappear when processing sparse matrices
  if(is.null(rownames(Xq))) {
    rownames(Xq) <- rownames(X)
  }
  #### Ensure sparse matrix ####
  Xq <- to_sparse(obj = Xq,
                  verbose = verbose)
  #### Compute specificity: step 1 ####
  Xq <- Matrix::t(Matrix::t(Xq) * (1 / Matrix::colSums(Xq)))
  #### Compute specificity: step 2 ####
  return(
    Xq / (apply(Xq, 1, sum) + 1e-12)
  )
}
