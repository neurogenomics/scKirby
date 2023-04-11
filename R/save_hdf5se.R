save_hdf5se <- function(obj,
                        save_dir,
                        quicksave_HDF5=TRUE,
                        overwrite=FALSE,
                        verbose=TRUE){
  if(!is_class(obj,"se")){
    stopper("obj must be a SummarizedExperiment object.")
  }

  if(!is.null(save_dir)){
    if(isTRUE(quicksave_HDF5) &&
       file.exists(file.path(save_dir,"assays.h5")) &&
       isFALSE(overwrite)){
      messager("+ Updating existing HDF5SummarizedExperiment:",save_dir,
               v=verbose)
      obj <- tryCatch(expr = {
        HDF5Array::quickResaveHDF5SummarizedExperiment(x=obj,
                                                       verbose=verbose)
      },
      error=function(e){
        messager("+ Saving HDF5SummarizedExperiment:",save_dir,v=verbose)
        hdf5_hard_save(obj=obj,
                       save_dir=save_dir,
                       overwrite=overwrite,
                       verbose=verbose)
      })
    } else {
      if( (!dir.exists(save_dir)) ||
          isTRUE(overwrite) ){
        messager("+ Saving HDF5SummarizedExperiment:",save_dir,v=verbose)
        obj <- hdf5_hard_save(obj=obj,
                              save_dir=save_dir,
                              overwrite=overwrite,
                              verbose=verbose)

      } else {
        messager("+ Returning existing SingleCellExperiment object.",v=verbose)
      }
    }
  }
  return(obj)
}
