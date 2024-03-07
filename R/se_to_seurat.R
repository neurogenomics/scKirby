#' Convert: \code{SummarizedExperiment} ==> \code{Seurat}
#'
#' @inheritParams converters
#' @export
#' @examples
#' obj <- example_obj("se")
#' obj2 <- se_to_seurat(obj)
se_to_seurat <- function(obj,
                         verbose=TRUE){
  messager_to()
  obj <- se_to_sce(obj = obj,
                   verbose = verbose)
  #### Ensure "counts" assay is present ####
  if(!"counts" %in% names(obj@assays)){
    names(obj@assays)[1] <- "counts"
  }
  #### Ensure "logcounts" assay is present ####
  if(!"logcounts" %in% names(obj@assays)){
    obj@assays@data$logcounts <- Seurat::LogNormalize(obj@assays@data$counts,
                                                      verbose = verbose)
  }
  Seurat::as.Seurat(obj)
}


