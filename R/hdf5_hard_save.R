hdf5_hard_save <- function(obj,
                           save_dir,
                           overwrite=TRUE,
                           verbose=TRUE){
  # DON'T create the HDF5 dir itself (will return an error about overwriting)
  # IMPORTANT!:
  ## set as.sparse=TRUE if you have the latest version of HDF5Array (1.8.11)
  if(packageVersion("HDF5Array")>=1.8){
    obj <- HDF5Array::saveHDF5SummarizedExperiment(x=obj,
                                                   dir=save_dir,
                                                   verbose=verbose,
                                                   as.sparse=TRUE,
                                                   replace=overwrite)
  }else{
    obj <- HDF5Array::saveHDF5SummarizedExperiment(x=obj,
                                                   dir=save_dir,
                                                   verbose=verbose,
                                                   replace=overwrite)
  }
  return(obj)
}
