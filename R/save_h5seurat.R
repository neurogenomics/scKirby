save_h5seurat <- function(obj,
                          save_path,
                          overwrite=TRUE,
                          verbose=TRUE){
  if(!is_class(obj,"seurat")){
    stopper("obj must be a Seurat object.")
  }
  messager("+ Saving h5Seurat:",save_path,v=verbose)
  obj_out <- SeuratDisk::as.h5Seurat(x = obj,
                                     filename = save_path,
                                     overwrite = overwrite)
  return(list(obj=obj_out,
              save_path=save_path))
}
