#' Extract reductions
#'
#' Extract dimensionality keys objects from any single-cell object class.
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
#' @returns A matrix of cell embeddings coordinates.
#'
#' @export
#' @importFrom methods is
#' @examples
#' obj <- example_obj("cds")
#' obsm <- get_obsm(obj)
get_obsm <- function(obj,
                     keys = NULL,
                     verbose = TRUE) {
  # devoptera::args2vars(get_obsm)

  #### Seurat: DimReduc ####
  if (methods::is(obj, "DimReduc")) {
    messager("Extracting obsm from DimReduc.", v = verbose)
    obsm <- dimreduc_to_list(r = obj)$obsm
  #### Seurat ####
  } else if (is_class(obj,"seurat")) {
    ## Seurat V1
    if(methods::is(obj,"seurat")){
      ### Need a way to figure which slot names are availabe a priori ####
      # messager("Extracting obsm from Seurat (V1).",v = verbose)
      # obsm <- list(PCA=list(embeddings=as.matrix(obj@pca.x),
      #                             loadings=as.matrix(obj@pca.rot)))
      obsm <- NULL;
    ## Seurat V2+
    } else {
      messager("Extracting obsm from Seurat.",v = verbose)
      all_keys <- Seurat::Reductions(obj)
      keys <- check_keys(all_keys = all_keys,
                              keys = keys)
      obsm <- lapply(stats::setNames(keys,keys),
                     function(nm){
                       dimreduc_to_list(r = obj@reductions[[nm]])$obsm
                       })
    }
  #### anndata ####
  } else if(is_class(obj,"anndata")){
    all_keys <- obj$obsm_keys()
    keys <- check_keys(
      all_keys = all_keys,
      keys = c(keys,
               if(is.null(keys)) {NULL} else {paste0("X_",keys)})
      )
    obsm <- lapply(stats::setNames(keys,keys),
           function(k){
             obj$obsm[[k]]
    })
  #### prcomp ####
  } else if (methods::is(obj, "prcomp")) {
    messager("Extracting obsm from prcomp.", v = verbose)
    obsm <- obj$x
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
    obsm <- if(methods::.hasSlot(obj,"reducedDimS") &&
               ncol(obj@reducedDimS)==ncol(obj) ){
      t(methods::slot(obj, "reducedDimS"))
    } else if (methods::.hasSlot(obj,"reducedDimA") &&
               ncol(obj@reducedDimA)==ncol(obj) ) {
      t(methods::slot(obj, "reducedDimA"))
    } else {
      NULL
    }
  #### list ####
  } else if (is_class(obj,"list")) {
    messager("Extracting obsm from list", v = verbose)
    if(is.character(obj$obsm)){
      obsm <- read_data(path = obj$obsm,
                        as_sparse = FALSE,
                        verbose = verbose)
    } else {
      obsm <- obj$obsm
    }
    #### Get keys ####
    if(!is.null(obsm)) {
      if(is.null(keys)) keys <- names(obsm)
      obsm <- obsm[keys]
    }
  } else if (is_class(obj,"matrix") ){
    messager("Interpretting obj as matrix with trait obsm.",
             v=verbose)
    obsm <- obj
  } else {
    stopper("No obsm could be identified.")
  }
  return(obsm)
}
