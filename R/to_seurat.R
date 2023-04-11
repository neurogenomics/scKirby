to_seurat <- function(obj,
                      workers=1,
                      verbose=TRUE){

  messager("Converting formats:",v=verbose)
  assign_cores(workers = workers)
  #### Check if class is supported ####
  check_supported(obj)
  #### CTD  ####
  if(is_class(obj,"ctd")){
    obj2 <- ctd_to_seurat(obj, verbose)
    return(obj2)
  }
  #### Matrices ####
  if(is_class(obj,"matrix")){
    obj2 <- matrix_to_seurat(obj, verbose)
    return(obj2)
  }
  #### Seurat ####
  if(is_class(obj,"seurat")){
    messager("+ Object already in Seurat format. Returning as-is.")
    return(obj)
  }
  #### anndata ####
  if(is_class(obj,"anndata")){
    obj2 <- anndata_to_seurat(obj,
                              verbose)
    return(obj2)
  }
  #### loom ####
  if(is_class(obj,"loom")){
    obj2 <- loom_to_seurat(obj, verbose)
    return(obj2)
  }
  #### SummarizedExperiment ####
  if(is_class(obj,"se")){
    obj2 <- sce_to_seurat(obj, verbose)
    return(obj2)
  }
  #### CellDataSet ####
  if(is_class(obj,"cds")){
    obj2 <- cds_to_seurat(obj, verbose)
    return(obj2)
  }
}
