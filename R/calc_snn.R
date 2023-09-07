calc_snn <- function(obj,
                     verbose=TRUE) {
    if (is.null(names(obj@graphs))) {
        if (!"pca" %in% names(obj@reductions)) {
            if (length(Seurat::VariableFeatures(obj)) == 0) {
                messager("No variable features detected. Computing",
                         v = verbose)
                obj <- Seurat::FindVariableFeatures(obj)
            }
            messager("No PCA detected. Computing", v = verbose)
            obj <- Seurat::NormalizeData(obj)
            obj <- Seurat::ScaleData(obj)
            obj <- Seurat::RunPCA(obj)
        }
        messager("No graphs detected. Computing.", v = verbose)
        obj <- Seurat::FindNeighbors(obj)
    }
    return(obj)
}
