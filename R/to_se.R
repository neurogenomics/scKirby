#' Convert: ==> \code{SummarizedExperiment}
#'
#' Convert any object to \code{SummarizedExperiment} or
#' \code{SingleCellExperiment} format.
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' l <- to_list(obj)
to_se <- function(obj,
                  as_sce=FALSE,
                  workers=1,
                  verbose=TRUE){
  messager("Converting formats:",v=verbose)
  assign_cores(workers = workers)
  #### Check if class is supported ####
  check_supported(obj)
  #### EWCE ####
  if(is_class(obj,"ewce")){
    obj2 <- ctd_to_se(obj,
                     as_sce = as_sce,
                     verbose = verbose)
    return(obj2)
  }
  #### Matrices ####
  if(is_class(obj,"matrix")){
    obj2 <- matrix_to_se(obj,
                        as_sce = as_sce,
                        verbose = verbose)
    return(obj2)
  }
  #### Seurat ####
  if(is_class(obj,"seurat")){
    obj2 <- seurat_to_se(obj,
                        as_sce = as_sce,
                        verbose = verbose)
    return(obj2)
  }
  #### AnnData ####
  if(is_class(obj,"anndata")){
    obj2 <- anndata_to_se(obj,
                          as_sce = as_sce,
                          verbose = verbose)
    return(obj2)
  }
  #### loom ####
  if(is_class(obj,"loom")){
    obj2 <- loom_to_se(obj, verbose)
    return(obj2)
  }
  #### SummarizedExperiment ####
  if(is_class(obj,"se")){
    messager("+ obj already in SummarizedExperiment format.",
             "Returning as-is.",v=verbose)
    obj2 <- check_se_rownames(obj, verbose = verbose)
    return(obj2)
  }
}
