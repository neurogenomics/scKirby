save_h5seurat <- function(object,
                          save_dir=tempdir(),
                          filename="scKirby_output",
                          verbose=T){
  filepath <- file.path(save_dir,filename)
  dir.create(dirname(filepath), showWarnings = F, recursive = T)
  messager("+ Seurat ==> h5Seurat",v=verbose)
  messager("+ Saving h5Seurat:",filepath,v=verbose)
  object_out <- SeuratDisk::as.h5Seurat(object,
                                        filename=filepath,
                                        overwrite=overwrite)
  return(list(object=object_out,
              filepath=filepath))
}
