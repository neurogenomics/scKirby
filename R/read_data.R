#' Read data
#'
#' Automatically infer the format of a single-cell dataset on disk and infer
#' how to read it into R.
#' @param path Path to saved single-cell data file,
#' or a directory containing the data files.
#' If \code{path} is already an R object, it will be returned directly.
#' @param filetype The type of file trying to be read in,
#'  made explicit by the user.
#'  "guess" (default) will simply infer the most likely file type.
#' @param custom_reader A user-supplied function to read in the data with.
#' @param verbose Print messages.
#' @inheritParams echoconda::activate_env
#' @returns An single-cell data object.
#'  The object class that the data gets imported as depends on file type.
#'
#' @export
#' @importFrom echoconda activate_env
#' @examples
#' library(Seurat)
#' obj <- example_obj("loom")
#' obj2 <- read_data(obj$filename)
read_data <- function(path,
                      filetype="guess",
                      custom_reader=NULL,
                      conda_env="r-reticulate",
                      verbose=TRUE,
                      ...){

  #### Use user-provided reader func ####
  if(!is.null(custom_reader)){
    messager("+ Reading in with custom_reader function.",v=verbose)
    obj <- custom_reader(path, ...)
    return(obj)
  }
  #### Return object directly ####
  if(!is.null(path) &&
     !is.na(path) &&
     !methods::is(path,"character")){
    messager("+ Returning object directly.",v=verbose)
    return(path)
  }
  if(!file.exists(path)){
    stopper("Cannot find file at path:",path)
  }
  #### Infer how to read in file ####
  messager("+ Reading from disk.",v=verbose)
  #### Generic RDS ####
  if(is_suffix(path,"rds") ||
     is_filetype(filetype,"rds")){
    messager("+ Reading in .rds file of unknown type.",v=verbose)
    readRDS(path,
            ...)
  }
  #### Generic RDS ####
  if(is_suffix(path,"rdata") ||
     is_filetype(filetype,"rdata")){
    messager("+ Reading in .rda file of unknown type.",v=verbose)
    load_rdata(path)
  }
  #### .mtx folder ####
  if((is_suffix(path,"matrix")  && dir.exists(path)) ||
     is_filetype(filetype,"matrix_dir")){
    read_matrix_dir(path = path,
                    verbose = verbose,
                    ...)
  }
  #### .mtx/.csv/.tsv matrix ####
  if((is_suffix(path,"matrix") && !dir.exists(path) )||
     is_filetype(filetype,"matrix") ||
     is_filetype(filetype,"data.table")){
    read_matrix(path = path,
                verbose = verbose,
                ...)
  }
  #### AnnData ####
  if(is_suffix(path,"anndata") ||
     is_filetype(filetype,"anndata")){
    read_anndata(path = path,
                 verbose = verbose,
                 conda_env = conda_env,
                 ...)
  }
  #### H5Seurat ####
  if(is_suffix(path,"h5seurat") ||
     is_filetype(filetype,"h5seurat")){
    read_h5seurat(path = path,
                  verbose = verbose,
                  ...)
  }
  #### HDF5Array SummarizedExperiment/SingleCellExperiment ####
  if((is_suffix(path,"h5") || is_filetype(filetype,"h5"))  &&
     dir.exists(obj) ){
    if(file.exists(file.path(path,"assays.h5")) &&
       file.exists(file.path(path,"se.rds")) ){
      read_h5(path = path,
              verbose = verbose,
              ...)
    }
  }
  #### Loom ####
  if(is_suffix(path,"loom") ||
     is_filetype(filetype,"loom")){
    read_loom(path = path,
              verbose = verbose,
              ...)
  }
}


