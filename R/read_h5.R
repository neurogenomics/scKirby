read_h5 <- function(path,
                    verbose,
                    ...){
  messager("+ HDF5Array format (.h5) detected.",
           "Importing as SingleCellExperiment.",v=verbose)
  HDF5Array::loadHDF5SummarizedExperiment(path,
                                          ...)
  # if( endsWith(tolower(obj), suffix=c(".h5")) & (!dir.exists(obj))){
  #     # When .h5 is a file, not a folder
  #     # obj <- "/Volumes/bms20/projects/neurogenomics-lab/live/GitRepos/model_celltype_conservation/raw_data/scRNAseq/LaManno2020/LaManno2020_sparse.h5"
  #     # rhdf5::h5closeAll()
  #     # h5_contents <- rhdf5::h5ls(obj, datasetinfo=F)
  #     # obj <- rhdf5::h5read(file =obj, name = "/")
  #     obj <- anndata::read_hdf(obj, key = "/")
  # }
}
