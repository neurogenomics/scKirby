#' Convert: \code{matrix} ==> \code{SummarizedExperiment}
#'
#' @inheritParams converters
#' @inheritParams to_se
#' @export
#' @examples
#' obj <- example_obj("matrix")
#' obj2 <- matrix_to_se(obj)
matrix_to_se <- function(obj,
                         as_sce=FALSE,
                         verbose=TRUE){
  messager_to()
  obj2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(raw = DelayedArray::DelayedArray(
      methods::as(as.matrix(obj), "sparseMatrix"))
      ),
  )
  if(isTRUE(as_sce)){
    obj2 <- se_to_sce(obj = obj2,
                      verbose = verbose)
  }
  obj2 <- check_se_rownames(obj2, verbose = verbose)
  return(obj2)
}
