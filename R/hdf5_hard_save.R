
hdf5_hard_save <- function(sce,
                           save_dir,
                           overwrite=T,
                           verbose=T){
  messager("+ Writing new HDF5...",v=verbose)
  # DON'T create the HDF5 dir itself (will return an error about overwriting)
  dir.create(dirname(save_dir), showWarnings = F, recursive = T)
  # IMPORTANT!: set as.sparse=T if you have the latest version of HDF5Array (1.8.11)
  pkg_ver <- packageVersion("HDF5Array")
  pkg_ver_split <- strsplit(as.character(pkg_ver),".", fixed = T)[[1]]
  pkg_V <- as.numeric(paste(pkg_ver_split[1], pkg_ver_split[2], sep="."))
  if(pkg_V>=1.8){
    sce <- HDF5Array::saveHDF5SummarizedExperiment(x=sce,
                                                   dir=save_dir,
                                                   verbose=verbose,
                                                   as.sparse=T,
                                                   replace=overwrite)
  }else{
    sce <- HDF5Array::saveHDF5SummarizedExperiment(x=sce,
                                                   dir=save_dir,
                                                   verbose=verbose,
                                                   replace=overwrite)
  }
  return(sce)
}
