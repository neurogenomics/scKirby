#' Example single-cell data object
#'
#' Create a small example object
#' of scRNA-seq data, in one of the following formats:
#' \itemize{
#' \item \link[SummarizedExperiment]{SummarizedExperiment}
#' \item \link[SeuratObject]{Seurat-class}
#' \item \link[SeuratObject]{Seurat-class}
#' \item \link[anndata]{AnnData}
#' }
#' ,
#'
#'
#' \emph{NOTE:} If you get the error
#'  \code{Error: ModuleNotFoundError: No module named 'anndata'},
#'  try running \link[anndata]{install_anndata} and
#'  then restarting your R session.
#'
#' Certain file types (loom, anndata) can't be stored a .rda
#'  because they need to be set up for on-disk reading.
#' @returns Single-cell object.
#'
#' @export
#' @importFrom utils data
#' @importClassesFrom SeuratObject Graph
#' @importFrom SeuratObject Graphs
#' @examples
#' se <- example_obj("se")
#' sce <- example_obj("sce")
#' seurat <- example_obj("seurat")
#'
#' library(Seurat) # necessary
#' loom <- example_obj("loom")
example_obj <- function(class=c("SummarizedExperiment",
                                "SingleCellExperiment",
                                "Seurat",
                                "h5Seurat",
                                "CellDataSet",
                                "matrix",
                                "loom",
                                "anndata"),
                        save_path = file.path(tempdir(),"example"),
                        verbose = TRUE){

    class <- tolower(class[[1]])

    #### matrix ####
    if(class %in% c("matrix","m")){
      obj <- SeuratObject::pbmc_small@assays$RNA@counts
    #### SummarizedExperiment ####
    } else if(class %in% c("summarizedexperiment","se")){
      obj <- seurat_to_se(SeuratObject::pbmc_small)
      obj <- sce_to_se(obj = obj, verbose = verbose)
    } else if(class %in% c("singellcellexperiment","sce")){
      obj <- seurat_to_se(SeuratObject::pbmc_small)
    #### SingleCellExperiment ####
    } else if(class %in% c("seurat","s")){
        obj <- SeuratObject::pbmc_small
    #### SeuratDisk ####
    } else if(class %in% c("h5seurat","h5s")){
      obj <- save_h5seurat(obj = SeuratObject::pbmc_small,
                           save_path =  paste0(save_path,".h5Seurat"),
                           verbose = verbose)
    #### CellDataSet ####
    } else if(class %in% c("celldataset","cds")){
      obj <- seurat_to_cds(obj = SeuratObject::pbmc_small)
    #### loom ####
    } else if(class %in% c("loom","l")){
      if(!is.null(save_path)){
        obj <- SeuratDisk::as.loom(x = SeuratObject::pbmc_small,
                                   overwrite = TRUE,
                                   filename = save_path)
      } else {
        obj <- SeuratDisk::as.loom(x = SeuratObject::pbmc_small)
      }
    #### anndata ####
    } else if(class %in% c("hdf5array","h5")){
      obj <- seurat_to_se(SeuratObject::pbmc_small)
      obj <- save_hdf5se(obj = obj,
                         save_dir = paste0(save_path,"_h5"),
                         verbose = verbose)
    }else if(class %in% c("anndata","ad")){
        obj <- example_anndata(save_path = save_path,
                               verbose = verbose)
    #### EWCE CelltypeDataset ####
    } else if(class %in% c("celltypedataset","celltypedata","ctd","ewce")){
      obj <- ewceData::ctd()
    #### ERROR ####
    } else {
        stp <- paste("class must be one of:",
                     paste("\n -",
                           shQuote(eval(formals(example_obj)$class)),
                           collapse = "")
                     )
        stop(stp)
    }
    return(obj)
}
