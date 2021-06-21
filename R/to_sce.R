
to_sce <- function(object,
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
    sce <- EWCEctd_to_sce(object, verbose)
    return(sce)
  }

  #### Check if class is supported ####
  if(!class(object)[1] %in% cdict$supported_classes){
    stop("Unsupported class detected: ",class(object),"\n\n",
         "Data object must be at least one of the following classes:\n\n",
         paste(supported_classes_print, collapse = ", "))
  }
  ### EWCE_list ####
  if(class(object)[1] == "EWCElist"){
    sce <- EWCElist_to_sce(object, verbose)
    return(sce)
  }

  #### Matrices ####
  if(class(object)[1] %in% cdict$matrix){
    sce <- matrix_to_sce(object, verbose)
    return(sce)
  }
  #### Seurat ####
  if(class(object)[1] %in% cdict$seurat){
    sce <- seurat_to_sce(object, verbose)
    return(sce)
  }
  #### AnnData ####
  if(class(object)[1] %in% cdict$anndata){
    sce <- anndata_to_sce(object, verbose)
    return(sce)
  }
  if(class(object)[1] %in% cdict$loom){
    sce <- loom_to_sce(object, verbose)
    return(sce)
  }
  #### SingleCellExperiment/SummarizedExperiment ####
  if(class(object)[1] %in% cdict$sce){
    messager("+ Object already in SingleCellExperiment format. Returning as-is.",v=verbose)
    sce <- check_sce_rownames(object, verbose = verbose)
    return(sce)
  }
}
