to_sparse <- function(obj,
                      allow_delayed_array = TRUE,
                      verbose = TRUE){

  if(is.null(rownames(obj)) &&
     is.character(obj[[1]])){
    obj <- as.matrix(obj[,-1], obj[[1]])
  }

  utils::getFromNamespace("to_sparse","orthogene")(
    gene_df2 = obj,
    allow_delayed_array = allow_delayed_array,
    verbose = verbose)
}
