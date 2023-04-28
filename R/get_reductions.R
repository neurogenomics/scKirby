#' Extract reductions
#'
#' Extract dimensionality keys objects from any single-cell object class.
#' @param obj Object of class:
#' \itemize{
#' \item{\code{DimReduc}}
#' \item{\code{Seurat}}
#' \item{\code{Seurat}}
#' }
#' @param keys The keys of reductions to extract from.
#' @param verbose Print messages.
#' @returns A named list with up to three elements:
#' \itemize{
#' \item{"embeddings": } Embeddings of samples (cells) in low-dimensional space.
#' \item{"loadings": } Loadings of features (genes) for each latent component.
#' \item{"loadings_projected": } Projected loadings of features (genes)
#'  for each latent component.
#' }
#'
#' @export
#' @importFrom methods is
#' @examples
#' obj <- example_obj("ad")
#' red <- get_reductions(obj)
get_reductions <- function(obj,
                           keys = NULL,
                           verbose = TRUE) {

  #### Seurat: DimReduc ####
  if (methods::is(obj, "DimReduc")) {
    messager("Extracting reductions from DimReduc.", v = verbose)
    reductions <- get_reductions_dimreduc(r = obj)
  #### Seurat ####
  } else if (is_class(obj,"seurat")) {
    ## Seurat V1
    if(methods::is(obj,"seurat")){
      messager("Extracting reductions from Seurat (V1).",v = verbose)
      reductions <- list(PCA=list(embeddings=obj@pca.x,
                                  loadings=obj@pca.rot))
    ## Seurat V2+
    } else {
      messager("Extracting reductions from Seurat.",v = verbose)
      all_keys <- Seurat::Reductions(obj)
      keys <- check_reduction_keys(all_keys = all_keys,
                                   keys = keys)
      reductions <- lapply(stats::setNames(keys,keys),
                           function(nm){
                             get_reductions_dimreduc(r = obj@reductions[[nm]])
                           })
    }
  #### anndata ####
  } else if(is_class(obj,"anndata")){
    all_keys <- obj$obsm_keys()
    keys <- check_reduction_keys(
      all_keys = all_keys,
      keys = c(keys,
               if(is.null(keys)) {NULL} else {paste0("X_",keys)})
      )
    reductions <- lapply(stats::setNames(keys,keys),
           function(k){
      list(embeddings = obj$obsm[[k]],
           loadings = obj$varm[[k]])
    })
  #### prcomp ####
  } else if (methods::is(obj, "prcomp")) {
    messager("Extracting reductions from prcomp.", v = verbose)
    reductions <- list(embeddings=obj$x,
                       loadings=obj$rotation)
  #### cds ####
  } else if (is_class(obj,"cds")){
    # nms <- grep("^reducedDim",methods::slotNames(obj),value = TRUE)
    key <- c(
      ## new coordinates of your cells in the reduced space.
      embeddings="reducedDimS",

      ## expression values for each cell (as a matrix) after
      ## whitening during dimensionality reduction.
      # reducedDimW="reducedDimW",

      ## Retrieves the the whitening matrix during
      ## independent component analysis.
      # reducedDimK="reducedDimK",

      ## The weights that transform the cells' coordinates in
      ## the reduced dimension space back to the full (whitened) space.
      loadings="reducedDimA")
    reductions <- lapply(key,function(nm){ methods::slot(obj, nm)})
  #### list ####
  } else if (is_class(obj,"list")) {
    messager("Extracting reductions from list", v = verbose)
    reductions <- obj$reductions
  } else if (is_class(obj,"matrix") ){
    messager("Interpretting obj as matrix with trait reductions.",
             v=verbose)
    reductions <- obj
  } else {
    stopper("No reductions could be identified.")
  }
  return(reductions)
}
