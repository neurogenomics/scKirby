#' Convert: \code{matrix} ==> \code{Seurat}
#'
#' @export
#' @examples
#' obj <- example_obj("list")
#' obj2 <- list_to_seurat(obj)
list_to_seurat <- function(obj,
                           verbose=TRUE){

  messager("+ list ==> Seurat",v=verbose)
  if(!is.null(obj$obs)){
    meta.data <- obj$obs
  } else {
    meta.data <- data.frame(cellid = colnames(obj$data),
                            row.names = colnames(obj$data)
                            )
  }
  aobj_list <- matrices_to_assayobjects(matrices = obj$data,
                                        var_features = obj$var_features,
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
  for(nm in names(obj$reductions)){
    r <- obj$reductions[[nm]]
    obj2[[nm]] <- Seurat::CreateDimReducObject(embeddings = r$embeddings,
                                               loadings = r$loadings,
                                               projected = r$loadings_projected,
                                               key = nm)
  }
  return(obj2)
}
