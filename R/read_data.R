
read_data <- function(obj,
                         filetype="guess",
                         custom_reader=NULL,
                         sce_save_dir=NULL,
                         overwrite=F,
                         verbose=T,
                         ...){
  if(!is.null(custom_reader)){
    messager("+ Reading in with custom_reader function...")
    object <- custom_reader(obj, ...)
    return(object)
  }

  if(class(obj)[1]=="character"){
    messager("+ Reading from disk...")
    #### Generic RDS ####
    if(endsWith(tolower(obj), suffix=".rds") | tolower(filetype)=="rds"){
      messager("+ Reading in .rds file of unknown type...")
      object <- readRDS(obj, ...)
      return(object)
    }
    #### Generic RDS ####
    if(any(endsWith(tolower(obj), suffix=c(".rda",".rdata"))) | tolower(filetype)=="rda"){
      messager("+ Reading in .rda file of unknown type...")
      object <- loadRData(obj)
      return(object)
    }
    #### .mtx folder ####
    if((endsWith(tolower(obj), suffix=".mtx") & dir.exists(obj)) | tolower(filetype)=="mtx_dir"){
      messager("+ Matrix folder format (.mtx) detected. Importing as Seurat object...")
      object <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(data.dir=obj, ...))
      return(object)
    }
    #### .mtx matrix ####
    if((endsWith(tolower(obj), suffix=".mtx") & (!dir.exists(obj)) ) | tolower(filetype)=="mtx"){
      messager("+ Expression matrix (.mtx) detected. Importing as sparse dgCMatrix...",v=verbose)
      object <- data.table::fread(obj, stringsAsFactors = F,
                                  data.table = F)
      object <- object %>% tibble::column_to_rownames(colnames(object)[1]) %>%
        as.matrix() %>% Matrix::Matrix(sparse=T)
      return(object)
    }
    #### .csv/.tsv matrix ####
    if((any(endsWith(tolower(obj), suffix=c(".csv",".tsv",".csv.gz",".tsv.gz")))) | tolower(filetype) %in% c("csv","tsv")){
      messager("+ Expression matrix (.csv|.tsv) detected. Importing as sparse dgCMatrix...",v=verbose)
      object <- data.table::fread(obj, stringsAsFactors = F,
                                  data.table = F)
      object <- object %>% tibble::column_to_rownames(colnames(object)[1]) %>%
        as.matrix() %>% Matrix::Matrix(sparse=T)
      return(object)
    }
    #### AnnData ####
    if((endsWith(tolower(obj), suffix=".h5ad")) | tolower(filetype)=="h5ad"){
      messager("+ AnnData format (.h5ad) detected. Importing as AnnData object...",v=verbose)
      #### anndata method
      # Anndata adds another dependency, but at least it works unlike
      object <- anndata::read_h5ad(filename = obj)

      #### sceasy method
      # object <- sceasy::convertFormat(obj, from="anndata", to="sce")

      #### Seurat method
      ## This is now deprecated, and for some reason doesn't offer back compatibility by calling to SeuratDisk...
      # object <- Seurat::ReadH5AD(obj, ...)

      ## SeuratDisk is currently broken, with no word from the developer...
      ## https://github.com/mojaveazure/seurat-disk/issues/41
      # object <- SeuratDisk::Convert(source = obj,
      #                               dest = sce_save_dir,
      #                               overwrite = overwrite,
      #                               verbose = verbose,
      #                               ...)
      return(object)
    }
    #### H5Seurat ####
    if((endsWith(tolower(obj), suffix=".h5seurat")) | tolower(filetype)=="h5seurat"){
      messager("+ h5Seurat format (.h5Seurat) detected. Importing as Seurat object...",v=verbose)
      object <- SeuratDisk::LoadH5Seurat(file = obj, ...)
      return(object)
    }
    #### Loom ####
    if((endsWith(tolower(obj), suffix=".loom")) | tolower(filetype)=="loom"){
      messager("+ Loom format (.loom) detected. Importing as SingleCellLoomExperiment object...",v=verbose)
      #### anndata method
      ## anndata::read_loom has difficulties identifying right loompy location.
      # anndata::read_loom(filename=obj, validate=F, ...)

      #### loomR method
      ## skip.validate must =F, or else you won't be able to extract the matrix
      # object <- loomR::connect(filename=obj, skip.validate = F)

      ### SeuratDisk method
      object <- SeuratDisk::LoadLoom(obj)

      #### sceasy method
      # rhdf5::h5disableFileLocking() ## Causes error otherwise
      # object <- sceasy::convertFormat(obj, from="loom", to="sce", ...)
      # outFile = gsub(".loom",".sce.rds",obj))
      return(object)
    }
    #### HDF5Array SummarizedExperiment/SingleCellExperiment ####
    if((any(endsWith(tolower(obj), suffix=c("h5","sce"))) |
        tolower(filetype) %in% c("HDF5Array","SummarizedExperiment","SingleCellExperiment","SingleCellLoomExperiment") ) &
       dir.exists(obj) ){
      if(file.exists(file.path(obj,"assays.h5")) & file.exists(file.path(obj,"se.rds")) ){
        messager("+ HDF5Array format (.h5) detected. Importing as SingleCellExperiment object...",v=verbose)
        object <- HDF5Array::loadHDF5SummarizedExperiment(obj, ...)
        return(object)
      }
    }

    # if( endsWith(tolower(obj), suffix=c(".h5")) & (!dir.exists(obj))){
    #     # When .h5 is a file, not a folder
    #     # obj <- "/Volumes/bms20/projects/neurogenomics-lab/live/GitRepos/model_celltype_conservation/raw_data/scRNAseq/LaManno2020/LaManno2020_sparse.h5"
    #     # rhdf5::h5closeAll()
    #     # h5_contents <- rhdf5::h5ls(obj, datasetinfo=F)
    #     # object <- rhdf5::h5read(file =obj, name = "/")
    #     object <- anndata::read_hdf(obj, key = "/")
    # }

  } else {
    messager("+ Returning object directly...",v=verbose)
    return(obj)
  }
}


