

#' Convert: \code{loom} ==> \code{SingleCellExperiment}
#'
#' @examples
#' library(scKirby)
#' sce <- loom_to_sce(object=example_loom())
#' @examples
loom_to_sce <- function(object,
                        verbose=T,
                        ...){
  messager("+ loom ==> SingleCellExperiment",v=verbose)
  # Import as a Seurat object first for convenience
  object <- SeuratDisk::LoadLoom(file = object$filename)
  # object <- loomR::connect(object$filename, mode = "r+")
  assay1 <- names(object@assays)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays      = list(raw = DelayedArray::DelayedArray(as(object@assays[[assay1]]@counts, "sparseMatrix"))),
    colData     = object@meta.data,
    rowData     = object@assays[[assay1]]@meta.features
  )
  sce <- check_sce_rownames(sce, verbose = verbose)
  return(sce)
}
