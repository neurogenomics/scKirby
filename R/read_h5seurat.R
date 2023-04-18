read_h5seurat <- function(path,
                          verbose=TRUE,
                          ...){
  messager("+ h5Seurat format (.h5Seurat) detected.",
           "Importing as Seurat.",v=verbose)
  SeuratDisk::LoadH5Seurat(file = path,
                           ...)
}
