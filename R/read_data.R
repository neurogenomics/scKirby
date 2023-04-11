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
#' @param verbose Print messages
#' @returns An single-cell data object.
#'  The object class that the data gets imported as depends on file type.
#'
#' @export
#' @examples
#' library(Seurat)
#' obj <- example_obj("loom")
#' obj2 <- read_data(obj$filename)
read_data <- function(path,
                      filetype="guess",
                      custom_reader=NULL,
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
    obj <- readRDS(path,
                   ...)
    return(obj)
  }
  #### Generic RDS ####
  if(is_suffix(path,"rdata") ||
     is_filetype(filetype,"rdata")){
    messager("+ Reading in .rda file of unknown type.")
    obj <- load_rdata(path)
    return(obj)
  }
  #### .mtx folder ####
  if((is_suffix(path,"matrix")  && dir.exists(path)) ||
     is_filetype(filetype,"matrix_dir")){
    messager("+ Matrix folder format (.mtx) detected.",
             "Importing as Seurat.")
    obj <- Seurat::CreateSeuratObject(
      counts = Seurat::Read10X(data.dir = path,
                               ...)
    )
    return(obj)
  }
  #### .mtx/.csv/.tsv matrix ####
  if((is_suffix(path,"matrix") && !dir.exists(path) )||
     is_filetype(filetype,"matrix") ||
     is_filetype(filetype,"data.table")){
    messager("+ Expression matrix (.mtx) detected.",
             "Importing as sparse matrix.",v=verbose)
    obj <- data.table::fread(path,
                             data.table = FALSE)
    obj2 <- as.matrix(obj[,-1], obj[[1]]) |>
      Matrix::Matrix(sparse=TRUE)
    return(obj2)
  }
  #### AnnData ####
  if(is_suffix(path,"anndata") ||
     is_filetype(filetype,"anndata")){
    messager("+ AnnData format (.h5ad) detected.",
             "Importing as AnnData.",v=verbose)
    #### anndata method
    # Anndata adds another dependency, but at least it works unlike
    obj <- anndata::read_h5ad(filename = path)

    #### sceasy method
    # obj <- sceasy::convertFormat(obj, from="anndata", to="sce")

    #### Seurat method
    ## This is now deprecated, and for some reason
    ## doesn't offer back compatibility by calling to SeuratDisk...
    # obj <- Seurat::ReadH5AD(obj, ...)

    ## SeuratDisk is currently broken, with no word from the developer...
    ## https://github.com/mojaveazure/seurat-disk/issues/41
    # obj <- SeuratDisk::Convert(source = obj,
    #                               dest = save_dir,
    #                               overwrite = overwrite,
    #                               verbose = verbose,
    #                               ...)
    return(obj)
  }
  #### H5Seurat ####
  if(is_suffix(path,"h5seurat") ||
     is_filetype(filetype,"h5seurat")){
    messager("+ h5Seurat format (.h5Seurat) detected.",
             "Importing as Seurat.",v=verbose)
    obj <- SeuratDisk::LoadH5Seurat(file = path,
                                    ...)
    return(obj)
  }
  #### HDF5Array SummarizedExperiment/SingleCellExperiment ####
  if((is_suffix(path,"h5") || is_filetype(filetype,"h5"))  &&
     dir.exists(obj) ){
    if(file.exists(file.path(path,"assays.h5")) &&
       file.exists(file.path(path,"se.rds")) ){
      messager("+ HDF5Array format (.h5) detected.",
               "Importing as SingleCellExperiment.",v=verbose)
      obj <- HDF5Array::loadHDF5SummarizedExperiment(path,
                                                     ...)
      return(obj)
    }
  }
  #### Loom ####
  if(is_suffix(path,"loom") ||
     is_filetype(filetype,"loom")){
    messager("+ Loom format (.loom) detected.",
             "Importing as SingleCellLoomExperiment.",v=verbose)
    #### anndata method
    ## anndata::read_loom has difficulties identifying right loompy location.
    # anndata::read_loom(filename=obj, validate=F, ...)

    #### loomR method
    ## skip.validate must =F, or else you won't be able to extract the matrix
    # obj <- loomR::connect(filename=obj, skip.validate = F)

    ### SeuratDisk method
    obj <- SeuratDisk::LoadLoom(path)

    #### sceasy method
    # rhdf5::h5disableFileLocking() ## Causes error otherwise
    # obj <- sceasy::convertFormat(obj, from="loom", to="sce", ...)
    # outFile = gsub(".loom",".sce.rds",obj))
    return(obj)
  }
  # if( endsWith(tolower(obj), suffix=c(".h5")) & (!dir.exists(obj))){
  #     # When .h5 is a file, not a folder
  #     # obj <- "/Volumes/bms20/projects/neurogenomics-lab/live/GitRepos/model_celltype_conservation/raw_data/scRNAseq/LaManno2020/LaManno2020_sparse.h5"
  #     # rhdf5::h5closeAll()
  #     # h5_contents <- rhdf5::h5ls(obj, datasetinfo=F)
  #     # obj <- rhdf5::h5read(file =obj, name = "/")
  #     obj <- anndata::read_hdf(obj, key = "/")
  # }
}


