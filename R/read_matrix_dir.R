read_matrix_dir <- function(path,
                            verbose,
                            ...){
  messager("+ Matrix directory format (.mtx) detected.",
           "Importing as Seurat.",v=verbose)
  Seurat::CreateSeuratObject(
    counts = Seurat::Read10X(data.dir = path,
                             ...)
  )
}
