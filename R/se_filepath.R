#' Get \code{SummarizedExperiment} file path
#'
#' Extract the file path from an
#' \link[SummarizedExperiment]{SummarizedExperiment} object. Only works if the
#' object has been stored on-disk.
#' @inheritParams converters
#' @returns File path, or NULL.
#'
#' @export
#' @examples
#' obj <- example_obj("h5")
#' f <- get_filepath_se(obj)
get_filepath_se <- function(obj){
  file_info <- DelayedArray::seed(SummarizedExperiment::assay(obj))
  if("filepath" %in% slotNames(file_info)){
    return(file_info@filepath)
  } else { return(NULL) }
}
