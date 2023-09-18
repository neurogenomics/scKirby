#' Convert: \code{MOFA} ==> \code{DimReducObject}
#'
#' @param obj A trained MOFA model.
#' @inheritParams get_obsm
#' @inheritParams converters
#' @inheritParams SeuratObject::CreateDimReducObject
#' @returns A named list of a \pkg{Seurat} \code{DimReducObject}.
#'
#' @export
mofa_to_dimreduc <- function(obj,
                             keys=NULL,
                             assay=NULL,
                             verbose=TRUE){

  messager_to_()
  if(is.null(keys)){
    keys <- tolower(c("mofa",paste0("mofa_",names(obj@dim_red))))
  }
  obsm <- get_obsm(obj = obj,
                   keys = keys,
                   verbose = verbose)
  varm <- get_varm(obj = obj,
                   keys = keys,
                   verbose = verbose)
  lapply(stats::setNames(keys,keys),
         function(k){
            messager("Converting",k,"to DimReducObject.",v=verbose)
            if (k=="mofa_umap") {
              dr <- Seurat::CreateDimReducObject(
                embeddings = obsm[[k]],
                key = k,
                assay = assay)
            }
            if (k=="mofa") {
              dr <- Seurat::CreateDimReducObject(
                embeddings = obsm[[k]],
                loadings = varm[[k]],
                key = k,
                assay = assay)
            }
            dr@misc[["variance_explained"]] <- obj@cache$variance_explained
            return(dr)
          })
}
