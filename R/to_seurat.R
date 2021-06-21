


to_seurat <- function(object,
                      save_dir=tempdir(),
                      verbose=T){
  messager("Converting formats:",v=verbose)
  core_allocation <- assign_cores(worker_cores = .90)
  cdict <- class_dictionary()
  #### EWCElist class ####
  if(is.EWCElist(object)){
    class(object) <- "EWCElist"
  }
  #### CTD class ####
  if(is.EWCEctd(object)){
    seurat <- EWCEctd_to_seurat(object, verbose)
    return(seurat)
  }
  #### Check if class is supported ####
  if(!class(object)[1] %in% cdict$supported_classes){
    stop("Unsupported class detected: ",class(object),"\n\n",
         "Data object must be at least one of the following classes:\n\n",
         paste(supported_classes_print, collapse = ", "))
  }
  ### EWCE_list ####
  if(class(object)[1] == "EWCElist"){
    seurat <- EWCElist_to_seurat(object, verbose)
    return(seurat)
  }
  #### Matrices ####
  if(class(object)[1] %in% cdict$matrix){
    seurat <- matrix_to_seurat(object, verbose)
    return(seurat)
  }
  #### Seurat ####
  if(class(object)[1] %in% cdict$seurat){
    printer("+ Object already in Seurat format. Returning as-is.")
    return(object)
  }
  #### AnnData ####
  if(class(object)[1] %in% cdict$anndata){
    seurat <- anndata_to_seurat(object, save_dir, verbose)
    return(seurat)
  }
  if(class(object)[1] %in% cdict$loom){
    seurat <- loom_to_seurat(object, verbose)
    return(seurat)
  }
  #### SingleCellExperiment/SummarizedExperiment ####
  if(class(object)[1] %in% cdict$sce){
    seurat <- sce_to_seurat(object, verbose)
    return(seurat)
  }
}
