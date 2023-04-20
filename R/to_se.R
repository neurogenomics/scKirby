#' Convert: ==> \code{SummarizedExperiment}
#'
#' Convert any object to \code{SummarizedExperiment} or
#' \code{SingleCellExperiment} format.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- to_se(obj)
to_se <- function(obj,
                  as_sce=FALSE,
                  verbose=TRUE){

  #### Check if class is supported ####
  check_supported(obj)
  #### EWCE ####
  if(is_class(obj,"ewce")){
    obj2 <- ctd_to_se(obj,
                     as_sce = as_sce,
                     verbose = verbose)
  #### Matrices ####
  } else if(is_class(obj,"matrix")){
    obj2 <- matrix_to_se(obj,
                        as_sce = as_sce,
                        verbose = verbose)
  #### Seurat ####
  } else if(is_class(obj,"seurat")){
    obj2 <- seurat_to_se(obj,
                        as_sce = as_sce,
                        verbose = verbose)
  #### AnnData ####
  } else if(is_class(obj,"anndata")){
    obj2 <- anndata_to_se(obj,
                          as_sce = as_sce,
                          verbose = verbose)
  #### loom ####
  } else if(is_class(obj,"loom")){
    obj2 <- loom_to_se(obj, verbose)
  #### SummarizedExperiment ####
  } else if(is_class(obj,"se")){
    messager("+ obj already in SummarizedExperiment format.",
             "Returning as-is.",v=verbose)
    obj2 <- check_se_rownames(obj, verbose = verbose)

  } else {
    l <- to_list(obj = obj,
                 verbose = verbose)
    obj2 <- list_to_se(obj = l,
                       as_sce = as_sce,
                       verbose = verbose)
  }
  #### Return ####
  return(obj2)
}
