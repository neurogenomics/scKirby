se_filepath <- function(obj){
  file_info <- DelayedArray::seed(SummarizedExperiment::assay(obj))
  if("filepath" %in% slotNames(file_info)){
    return(file_info@filepath)
  } else { return(NULL) }
}
