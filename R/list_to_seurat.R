#' Convert: \code{list} ==> \code{Seurat}
#'
#' @inheritParams converters
#' @export
#' @examples
#' obj <- example_obj("list")
#' obj2 <- list_to_seurat(obj)
list_to_seurat <- function(obj,
                           verbose=TRUE){

  messager_to()
  obs <- get_obs(obj = obj,
                 verbose = verbose)
  X_list <- get_x(obj)
  aobj_list <- matrices_to_assayobjects(X_list = X_list,
                                        var_features = obj$varm,
                                        verbose = verbose)
  obj2 <- SeuratObject::CreateSeuratObject(
    counts = aobj_list[[1]],
    assay = names(aobj_list)[[1]],
    meta.data = obs)
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
    messager("Creating DimReducObject:",shQuote(nm),v=verbose)
    obj2[[nm]] <- SeuratObject::CreateDimReducObject(
      embeddings = obsm[[nm]],
      loadings = if(is.null(varm[[nm]])) new(Class = "matrix") else varm[[nm]],
      assay = names(aobj_list)[[1]],
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
