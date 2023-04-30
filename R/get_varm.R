#' Extract variable loadings
#'
#' Extract variable feature loadings from any single-cell object class.
#' @param obj Object of class:
#' \itemize{
#' \item{\code{DimReduc}}
#' \item{\code{Seurat}}
#' \item{\code{AnnData}}
#' \item{\code{prcomp}}
#' \item{\code{CellDataSet}}
#' \item{\code{list}}
#' \item{\code{matrix/data.frame}}
#' }
#' @param keys The keys of reductions to extract from.
#' @param verbose Print messages.
#' @returns A matrix of feature embeddings coordinates.
#'
#' @export
#' @importFrom methods is
#' @examples
#' obj <- example_obj("ad")
#' varm <- get_varm(obj)
get_varm <- function(obj,
                     keys = NULL,
                     verbose = TRUE) {

  #### Seurat: DimReduc ####
  if (methods::is(obj, "DimReduc")) {
    messager("Extracting varm from DimReduc.", v = verbose)
    varm <- dimreduc_to_list(r = obj)$varm
  #### Seurat ####
  } else if (is_class(obj,"seurat")) {
    ## Seurat V1
    if(methods::is(obj,"seurat")){
      ### Need a way to figure which slot names are availabe a priori ####
      # messager("Extracting varm from Seurat (V1).",v = verbose)
      # varm <- list(PCA=list(embeddings=as.matrix(obj@pca.x),
      #                             loadings=as.matrix(obj@pca.rot)))
      varm <- NULL;
    ## Seurat V2+
    } else {
      messager("Extracting varm from Seurat.",v = verbose)
      all_keys <- Seurat::Reductions(obj)
      keys <- check_keys(all_keys = all_keys,
                         keys = keys)
      varm <- lapply(stats::setNames(keys,keys),
                     function(nm){
                       dimreduc_to_list(r = obj@reductions[[nm]])$varm
                     })
    }
  #### anndata ####
  } else if(is_class(obj,"anndata")){
    all_keys <- obj$varm_keys()
    keys <- check_keys(
      all_keys = all_keys,
      keys = c(keys,
               if(is.null(keys)) {NULL} else {paste0("X_",keys)})
    )
    varm <- lapply(stats::setNames(keys,keys),
                   function(k){
                     obj$varm[[k]]
                   })
  #### prcomp ####
  } else if (methods::is(obj, "prcomp")) {
    messager("Extracting varm from prcomp.", v = verbose)
    varm <- obj$rotation
  #### cds ####
  } else if (is_class(obj,"cds")){
    # key <- c(
    #   ## new coordinates of your cells in the reduced space.
    #   embeddings="reducedDimS",
    #
    #   ## expression values for each cell (as a matrix) after
    #   ## whitening during dimensionality reduction.
    #   # reducedDimW="reducedDimW",
    #
    #   ## Retrieves the the whitening matrix during
    #   ## independent component analysis.
    #   # reducedDimK="reducedDimK",
    #
    #   ## The weights that transform the cells' coordinates in
    #   ## the reduced dimension space back to the full (whitened) space.
    #   loadings="reducedDimA")
    varm <- if (methods::.hasSlot(obj,"reducedDimA") &&
               ncol(obj@reducedDimA)==nrow(obj) ) {
      t(methods::slot(obj, "reducedDimA"))
    } else {
      NULL
    }
  #### list ####
  } else if (is_class(obj,"list")) {
    messager("Extracting varm from list", v = verbose)
    varm <- obj$varm
  } else if (is_class(obj,"matrix") ){
    messager("Interpretting obj as matrix with trait varm.",
             v=verbose)
    varm <- obj
  } else {
    stopper("No varm could be identified.")
  }
  return(varm)
}
