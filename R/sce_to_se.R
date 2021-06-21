sce_to_se <- function(object,
                      verbose=T){
  messager("+ SingleCellExperiment ==> SummarizedExperiment",v=verbose)
  object_out <- SummarizedExperiment::SummarizedExperiment(object)
  return(object_out)
}
