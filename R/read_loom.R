read_loom <- function(path,
                      method = c("anndata","seuratdisk","scopeloomr"),
                      verbose = TRUE,
                      ...){

  method <- tolower(method)[1]

  #### anndata method ####
  if(method=="anndata"){
    activate_conda(verbose=verbose)
    # anndata::read_loom has difficulties identifying right loompy location.
    anndata::read_loom(filename=path,
                       validate=FALSE,
                       ...)

  ### SeuratDisk method #####
  } else if(method=="scopeloomr"){

    SCopeLoomR::open_loom(path)
  }else if(method=="seuratdisk"){
    messager("+ Loom format (.loom) detected.",
             "Importing as SingleCellLoomExperiment.",v=verbose)
    SeuratDisk::LoadLoom(file = path,
                         verbose = verbose,
                         ...)
  }
  #### loomR method
  ## loomR has since been superseded by SeuratDisk
  ## skip.validate must =F, or else you won't be able to extract the matrix
  # obj <- loomR::connect(filename=obj, skip.validate = F)

  #### sceasy method
  # rhdf5::h5disableFileLocking() ## Causes error otherwise
  # obj <- sceasy::convertFormat(obj, from="loom", to="sce", ...)
  # outFile = gsub(".loom",".sce.rds",obj))
}
