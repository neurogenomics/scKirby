
sce_filepath <- function(sce,
                         sce_save_dir=NULL){
  file_info <- DelayedArray::seed(SummarizedExperiment::assay(sce))
  if("filepath" %in% slotNames(file_info)){
    return(file_info@filepath)
  } else { return(NULL) }
}
