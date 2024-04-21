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
#' @inheritParams get_n_elements
#' @returns A named list of matrices containing the feature loadings.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' varm <- get_varm(obj)
get_varm <- function(obj,
                     keys = NULL,
                     n = NULL,
                     verbose = TRUE) {

  #### Matrix ####
  if (is_class(obj,"matrix") ){
    messager("Interpretting obj as matrix with trait varm.",
             v=verbose)
    varm <- list(varm=obj)
  #### MOFA2: model ####
  } else if(methods::is(obj,"MOFA")){
    messager("Extracting varm from MOFA", v = verbose)
    varm <- get_varm_mofa(obj = obj,
                          verbose = verbose)
  #### Seurat: DimReduc ####
  } else if (methods::is(obj, "DimReduc")) {
    messager("Extracting varm from DimReduc.", v = verbose)
    varm <- dimreduc_to_list(obj = obj)$varm
  #### Seurat ####
  } else if (is_class(obj,"seurat")) {
    ## Seurat V1
    if(methods::is(obj,"seurat")){
      ### Need a way to figure which slot names are available a priori ####
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
                       dimreduc_to_list(obj = obj@reductions[[nm]])$varm
                     })
    }
  #### AssayData ####
  } else if (methods::is(obj,"AssayData")) {
      varm <- if ("pca_models" %in% names(obj)) {
        list(pca=obj$pca_models[[1]]$rotation)
      } else if ("v" %in% names(obj)) {
        list(v=obj["v"])
      } else if ("V" %in% names(obj)) {
        list(V=obj$V)
      } else {
        messager("No varm could be identified.","Returning NULL.", v=verbose)
        NULL
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
    varm <- list(pca=obj$rotation)
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
      list(reducedDimA=t(methods::slot(obj, "reducedDimA")))
    } else {
      return(NULL)
    }
  #### list ####
  } else if (is_class(obj,"matrix_list")){
    messager("Interpretting obj as matrix list of varm.",v=verbose)
    return(obj)
  } else if (is_class(obj,"list")) {
    messager("Extracting varm from list", v = verbose)
    if(is.character(obj$varm)){
      varm <- read_data(path = obj$varm,
                        as_sparse = FALSE,
                        verbose = verbose)
    } else {
      varm <- obj[["varm"]]
    }
  } else {
    messager("No varm could be identified.","Returning NULL.", v=verbose)
    return(NULL)
  }
  if(is.null(NULL)) return(NULL)
  #### Return as a named list (1 per assay), unless there's only 1 assay ####
  varm <- get_n_elements(l = varm,
                         n = n,
                         verbose = verbose)
  return(varm)
}
