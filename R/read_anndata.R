read_anndata <- function(path,
                         method = c("anndata","zellkonverter"),
                         verbose=TRUE,
                         ...){

  method <- tolower(method)[1]
  messager("+ AnnData format (.h5ad) detected.",
           "Importing as AnnData.",v=verbose)
  #### Expand path for python ####
  if(!is.null(path)){
    path <- path.expand(path)
  }
  if(!file.exists(path)){
    stopper("File does not exist:",path)
  }
  #### Method: anndata ####
  if(method=="anndata"){
    #### anndata method
    # Anndata adds another dependency, but at least it works unlike
    activate_conda(verbose=verbose)
    obj <- anndata::read_h5ad(filename = path,
                              ...)
  #### Method: zellkonverter ####
  } else {
    requireNamespace("zellkonverter")
    #### Has to import as SCE first,
    ## but has more robust reading/conversion handling
    ## than other methods. Therefore mostly useful for fringe cases.
    obj <- zellkonverter::readH5AD(file = path,
                                   verbose = verbose,
                                   ...)
    obj <- zellkonverter::SCE2AnnData(sce = obj,
                                      verbose = verbose,
                                      ...)
  }
  return(obj)

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
  #
}
