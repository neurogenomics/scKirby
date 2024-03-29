#' Process a \code{Seurat} object
#'
#' Run a standardised \link[SeuratObject]{Seurat} pipeline on
#' a \link[SeuratObject]{Seurat} object
#' or raw \code{counts} and \code{meta.data}.\cr
#' Automatically performs
#' \describe{
#' \item{\code{FindVariableFeatures}}{Variable feature selection}
#' \item{\code{NormalizeData}}{Data normalization}
#' \item{\code{RunPCA}}{PCA}
#' \item{\code{RunUMAP}}{UMAP}
#' \item{\code{FindNeighbors}}{K-nearest neighbors}
#' \item{\code{FindClusters}}{Clustering}
#' }
#' @param cluster_reduction Recompute neighbors graph based on UMAP
#' to get clusters that best reflect UMAP space.
#' For this same reason, only cluster in two dimensions,
#' because this is the view we most often use.
#' That said, this may reduce the generalisability of these clusters/graph.
#' @param log_norm Log-normalise the data with \link[Seurat]{LogNormalize}.
#' @param max_mem Max limit on memory usage, to be passed to environmental
#'  variable \code{future.globals.maxSize}.
#' @param seed Random seed to set.
#' @param add_specificity Add a new assay called "specificity" by deriving
#' specificity scores from the input data.
#' @param default_assay Default assay to set with
#'  \link[SeuratObject]{DefaultAssay}.
#' @param workers Number of workers (threads) to parallelise processes across
#' using \pkg{future}.
#' @inheritParams converters
#' @inheritParams SeuratObject::CreateSeuratObject
#' @inheritParams Seurat::FindVariableFeatures
#' @inheritParams Seurat::ScaleData
#' @inheritParams Seurat::RunUMAP
#' @returns A preprocessed \link[SeuratObject]{Seurat} object.
#'
#' @export
#' @importFrom Seurat CreateSeuratObject GetAssayData CreateAssayObject
#' DefaultAssay FindVariableFeatures NormalizeData LogNormalize ScaleData
#' RunPCA FindNeighbors RunUMAP FindClusters
#' @importFrom future plan
#' @examples
#' obj <- example_obj("list")
#' obj2 <- process_seurat(obj = obj)
process_seurat <- function(obj = NULL,
                           meta.data = NULL,
                           nfeatures = 2000,
                           vars.to.regress = NULL,
                           dims = seq(50),
                           add_specificity = FALSE,
                           default_assay = NULL,
                           n.components = 2L,
                           log_norm = FALSE,
                           cluster_reduction = "umap",
                           workers = 1,
                           max_mem = 8000*1024^2,
                           seed = 2020,
                           verbose = TRUE){
  requireNamespace("future")
  set.seed(seed)

  future::plan(strategy = "multicore", workers = workers)
  options(future.globals.maxSize = max_mem)
  #### Ensure in Seurat format ####
  obj2 <- to_seurat(obj = obj,
                    verbose = verbose)
  #### Create new assay using EWCE-style specificity ####
  if(isTRUE(add_specificity)){
    Xs <- calc_specificity(Seurat::GetAssayData(obj2))
    obj2[['specificity']] <- Seurat::CreateAssayObject(counts = Xs)
    remove(Xs)
  }
  #### Set assay ####
  if(!is.null(default_assay)) {
    Seurat::DefaultAssay(obj2) <- default_assay
  }
  #### Select variable features ####
  if(is.null(nfeatures)) nfeatures <- nrow(obj2)
  obj2 <- Seurat::FindVariableFeatures(obj2,
                                       nfeatures = nfeatures)
  obj2 <- Seurat::NormalizeData(obj2)
  # if(isTRUE(log_norm)) logged <- Seurat::LogNormalize(obj2)
  obj2 <- Seurat::ScaleData(obj2,
                            vars.to.regress = vars.to.regress)

  #### Dimensionality reduction ####
  obj2 <- Seurat::RunPCA(obj2)
  obj2 <- Seurat::FindNeighbors(obj2)
  obj2 <- Seurat::RunUMAP(obj2,
                          dims = dims,
                          n.components = n.components,
                          return.model = TRUE)
  #### Clustering ####
  if(!isFALSE(cluster_reduction)){
    obj2 <- Seurat::FindNeighbors(obj2,
                                  reduction = cluster_reduction,
                                  dims = seq_len(n.components))
    obj2 <- Seurat::FindClusters(obj2,
                                 reduction = cluster_reduction)
  } else {
    obj2 <- Seurat::FindClusters(obj2)
  }
  #### Return ####
  return(obj2)
}
