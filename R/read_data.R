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
                      verbose=TRUE,
                      ...){
  # devoptera::args2vars(read_data)

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
  #### Infer how to read in file ####
  if(!file.exists(path)){
    stopper("Cannot find file at path:",path)
  }
  messager("+ Reading from disk.",v=verbose)
  #### Generic RDS ####
  if(is_suffix(path,"rds") ||
     is_filetype(filetype,"rds")){
    messager("+ Reading in .rds file of unknown type.",v=verbose)
    obj <- readRDS(path,
                   ...)
  #### Generic RDS ####
  } else if(is_suffix(path,"rdata") ||
     is_filetype(filetype,"rdata")){
    messager("+ Reading in .rda file of unknown type.",v=verbose)
    obj <- load_rdata(path)
  #### .mtx folder ####
  } else if((is_suffix(path,"matrix")  && dir.exists(path)) ||
     is_filetype(filetype,"matrix_dir")){
    obj <- read_matrix_dir(path = path,
                           verbose = verbose,
                           ...)
  #### .mtx/.csv/.tsv matrix ####
  } else if((is_suffix(path,"matrix") && !dir.exists(path) )||
            is_suffix(path,"data.table") ||
            is_filetype(filetype,"matrix") ||
            is_filetype(filetype,"data.table")){
    obj <- read_matrix(path = path,
                       verbose = verbose,
                       ...)
  #### AnnData ####
  } else if(is_suffix(path,"anndata") ||
     is_filetype(filetype,"anndata")){
    obj <- read_anndata(path = path,
                        verbose = verbose,
                        ...)
  #### H5Seurat ####
  } else if(is_suffix(path,"h5seurat") ||
     is_filetype(filetype,"h5seurat")){
    obj <- read_h5seurat(path = path,
                         verbose = verbose,
                         ...)
  #### HDF5Array SummarizedExperiment/SingleCellExperiment ####
  } else if((is_suffix(path,"h5") || is_filetype(filetype,"h5"))  &&
     dir.exists(obj) ){
    if(file.exists(file.path(path,"assays.h5")) &&
       file.exists(file.path(path,"se.rds")) ){
      obj <- read_h5(path = path,
                     verbose = verbose,
                     ...)
    }
  #### Loom ####
  } else if(is_suffix(path,"loom") ||
     is_filetype(filetype,"loom")){
    obj <- read_loom(path = path,
                     verbose = verbose,
                     ...)
  } else {
    stopper("Unable to read in file.")
  }
  return(obj)
}


