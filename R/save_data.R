#' Save data
#'
#' Save a variety of single-cell data formats.
#' @inheritParams ingest_data
#' @export
#' @examples
#' obj <- example_obj("cds")
#' filepath <- save_data(obj)
save_data <- function(obj,
                      filetype = c("h5","h5seurat","h5ad","rda","rds","anndata"),
                      save_path = file.path(tempdir(),
                                          paste("scKirby_output",

                                                filetype,sep=".")),
                      quicksave_HDF5 = TRUE,
                      overwrite = TRUE,
                      verbose = TRUE){

  #### Setup dir ####
  dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
  #### hdf5se ####
  if(is_filetype(filetype,"h5")){
    save_path <- save_hdf5se(obj=obj,
                             save_dir=save_path,
                             quicksave_HDF5=quicksave_HDF5,
                             overwrite=overwrite,
                             verbose=verbose)
  #### h5seurat ####
  } else if(is_filetype(filetype,"h5seurat")){
    object_filepath <- save_h5seurat(obj=obj,
                                     save_path=save_path,
                                     verbose=verbose)
    save_path <- object_filepath$save_path
  #### h5ad ####
  } else if(is_filetype(filetype,"anndata")){
    save_anndata(obj = obj,
                 save_path = save_path,
                 verbose = verbose)
  #### rda ####
  } else if(is_filetype(filetype,"rda")){
    messager("+ Saving RData:",save_path,v=verbose)
    save(obj,file = save_path)
  #### rds ####
  } else {
    messager("+ Saving RDS:",save_path,v=verbose)
    saveRDS(obj,save_path)
  }
  return(save_path)
}
