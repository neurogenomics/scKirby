

#' Convert: \code{Seurat} ==> \code{SingleCellExperiment}
#'
#' @examples
#' library(scKirby)
#' data("example_seurat")
#' sce <- seurat_to_sce(example_seurat)
seurat_to_sce <- function(object,
                          verbose=T){
  messager("+ Seurat ==> SingleCellExperiment",v=verbose)
  # object <- Seurat::pbmc_small ## example
  ### Seurat::as.SingleCellExperiment() is currently broken
  # sce <- Seurat::as.SingleCellExperiment(object)

  ### sceasy seurat --> sce is consequently also broken
  # sce <- sceasy::convertFormat(object, from = "seurat", to="sce")

  ### Must convert manually instead
  assay1 <- names(object@assays)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays      = list(raw = DelayedArray::DelayedArray(as(object@assays[[assay1]]@counts, "sparseMatrix"))),
    colData     = object@meta.data,
    rowData     = object@assays[[assay1]]@meta.features
  )
  sce <- check_sce_rownames(sce, verbose = verbose)
  return(sce)
}
