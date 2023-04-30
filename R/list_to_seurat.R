#' Convert: \code{list} ==> \code{Seurat}
#'
#' @export
#' @examples
#' obj <- example_obj("list")
#' obj2 <- list_to_seurat(obj)
list_to_seurat <- function(obj,
                           verbose=TRUE){

  #### Activate conda env with anndata installed ####
  activate_conda(verbose=verbose)
  messager("+ list ==> Seurat",v=verbose)
  if(!is.null(obj$obs)){
    meta.data <- obj$obs
  } else {
    meta.data <- data.frame(cellid = colnames(obj$data),
                            row.names = colnames(obj$data)
                            )
  }
  aobj_list <- matrices_to_assayobjects(matrices = obj$data,
                                        var_features = obj$varm,
                                        verbose = verbose)
  obj2 <- SeuratObject::CreateSeuratObject(
    counts = aobj_list[[1]],
    assay = names(aobj_list)[[1]],
    meta.data = obj$obs)
  #### Add extra assays ####
  if(length(aobj_list)>1){
    for(nm in names(aobj_list)[-1])
      obj2[[nm]] <- aobj_list[[nm]]
  }
  #### Add reductions ####
  obsm <- get_obsm(obj = obj,
                   verbose = verbose)
  varm <- get_varm(obj = obj,
                   verbose = verbose)
  for(nm in names(obj$obsm)){
    obj2[[nm]] <- Seurat::CreateDimReducObject(
      embeddings = obsm[[nm]],
      loadings = varm[[nm]],
      key = nm)
  }
  #### Add variable feature ####
  if("var_features" %in% names(obj$uns)){
    for(a in names(obj2@assays)){
      obj2@assays[[a]]@var.features <- obj$uns$var_features[[a]]
    }
  }
  return(obj2)
}
