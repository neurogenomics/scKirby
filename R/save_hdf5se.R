save_hdf5se <- function(sce,
                        save_dir,
                        quicksave_HDF5=T,
                        overwrite=F,
                        verbose=T){
  if(!is.null(save_dir)){
    if(quicksave_HDF5 & file.exists(file.path(save_dir,"assays.h5")) & (overwrite==F)){
      messager("+ Updating existing HDF5SummarizedExperiment:",save_dir,v=verbose)
      sce <- tryCatch(expr = {
        HDF5Array::quickResaveHDF5SummarizedExperiment(x=sce, verbose=verbose)
      },
      error=function(e){
        messager("+ Saving HDF5SummarizedExperiment:",save_dir,v=verbose)
        hdf5_hard_save(sce=sce,
                       save_dir=save_dir,
                       overwrite=overwrite,
                       verbose=verbose)
      })
    } else {
      if( (!dir.exists(save_dir)) | (overwrite) ){
        messager("+ Saving HDF5SummarizedExperiment:",save_dir,v=verbose)
        sce <- hdf5_hard_save(sce=sce,
                              save_dir=save_dir,
                              overwrite=overwrite,
                              verbose=verbose)

      } else {
        messager("+ Returning existing SingleCellExperiment object.",v=verbose)
      }
    }
  }
  return(sce)
}
