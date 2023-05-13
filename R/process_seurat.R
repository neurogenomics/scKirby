#' Process a \link[SeuratObject]{Seurat} object
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
#' @returns A preprocessed \link[Seurat]{Seurat} object.
#'
#' @export
#' @import Seurat
#' @import future
#' @examples
#' obj <- example_obj("matrix")
#' obj2 <- seurat_pipeline(obj = obj)
process_seurat <- function(obj = NULL,
                           meta.data = NULL,
                           nfeatures = 2000,
                           vars.to.regress = NULL,
                           dims = seq_len(50),
                           add_specificity = FALSE,
                           assay_name = "RNA",
                           default_assay = NULL,
                           n.components = 2L,
                           log_norm = FALSE,
                           cluster_reduction = "umap",
                           workers = 1,
                           max_mem = 8000*1024^2,
                           seed = 2020){
  requireNamespace("future")
  set.seed(seed)

  future::plan(strategy = "multicore", workers = workers)
  options(future.globals.maxSize = max_mem)

  if(!methods::is(obj,"Seurat")){
    obj2 <- Seurat::CreateSeuratObject(counts = obj,
                                       meta.data = meta.data,
                                       assay = assay_name)
  } else {
    obj2 <- obj
  }
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
  if(isTRUE(log_norm)) logged <- Seurat::LogNormalize(obj2)
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
