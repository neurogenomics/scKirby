


save_data <- function(object,
                      output_type,
                      save_dir=tempdir(),
                      filename="scKirby_output",
                      quicksave_HDF5=T,
                      overwrite=T,
                      verbose=T){
  cdict <- class_dictionary()
  if(tolower(output_type)  %in% tolower(cdict$sce)){
    filepath <- save_sce(object=object,
                         save_dir=save_dir,
                         filename=filename,
                         verbose=verbose)
  }
  if(tolower(output_type)  %in% tolower(cdict$hdf5se)){
    filepath <- save_hdf5se(sce=object,
                            save_dir=save_dir,
                            quicksave_HDF5=quicksave_HDF5,
                            overwrite=overwrite,
                            verbose=verbose)
  }
  if(tolower(output_type)  %in% tolower(cdict$seurat)){
    filepath <- save_seurat(object=object,
                            save_dir=save_dir,
                            filename=filename,
                            verbose=verbose)
  }
  if(tolower(output_type) %in% tolower(cdict$h5seurat)){
    object_filepath <- save_h5seurat(object=object,
                                     save_dir=save_dir,
                                     filename=filename,
                                     verbose=verbose)
    filepath <- object_filepath$filepath
    }
  return(filepath)
}
