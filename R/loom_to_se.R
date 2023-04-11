#' Convert: \code{loom} ==> \code{SummarizedExperiment}
#'
#' @export
#' @examples
#' library(Seurat)
#' obj <- example_obj("loom")
#' obj2 <- loom_to_se(obj)
loom_to_se <- function(obj,
                       as_sce=FALSE,
                       verbose=TRUE,
                        ...){
  messager("+ loom ==> SummarizedExperiment",v=verbose)
  #### Import as a Seurat obj first for convenience ####
  obj <- SeuratDisk::LoadLoom(file = obj$filename)
  #### Then convert to se/sce ####
  obj2 <- seurat_to_se(obj,
                       as_sce = as_sce,
                       verbose = verbose)
  return(obj2)
}
