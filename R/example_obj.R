#' Example single-cell data object
#'
#' Create a small example object
#' of scRNA-seq data.
#'
#' \emph{NOTE:} If you get the error
#'  \code{Error: ModuleNotFoundError: No module named 'anndata'},
#'  try running \link[anndata]{install_anndata} and
#'  then restarting your R session.
#'  Also, loom files require you to run \code{library(Seurat)} beforehand.
#'
#' Certain file types (loom, anndata) can't be stored as \emph{.rda} files
#'  because they need to be set up for on-disk reading.
#' @param class Object class to return.
#' @param save_path Path to save on-disk object formats.
#' @inheritParams converters
#' @returns A single-cell object.
#'
#' @export
#' @examples
#' se <- example_obj("se")
#' sce <- example_obj("sce")
#' seurat <- example_obj("seurat")
#'
#' library(Seurat) # necessary
#' loom <- example_obj("loom")
example_obj <- function(class=c("Seurat",
                                "h5Seurat",
                                "SummarizedExperiment",
                                "SingleCellExperiment",
                                "CellDataSet",
                                "matrix",
                                "data.table",
                                "data.frame",
                                "list",
                                "loom",
                                "LoomExperiment",
                                "SingleCellLoomExperiment",
                                "hdf5se",
                                "anndata",
                                "EWCE"),
                        save_path = file.path(tempdir(),"example"),
                        verbose = TRUE,
                        ...){
  # devoptera::args2vars(example_obj)

  class <- tolower(class[[1]])
  base_obj <- SeuratObject::pbmc_small
  #### matrix ####
  if(class %in% c("matrix","m")){
    obj <- get_x(obj = base_obj,
                 n = 1,
                 verbose = verbose,
                 ...)
  #### data.frame ####
  } else  if(class %in% c("matrix_list","ml")){
    obj <- get_x(obj = base_obj,
                 verbose = verbose,
                 ...)
  } else if(class %in% c("data.frame","df")){
    obj <- to_dataframe(obj = obj,
                        verbose = verbose,
                        ...)
  #### data.table ####
  } else if(class %in% c("data.table","dt")){
    obj <- to_datatable(obj = obj,
                        verbose = verbose,
                        ...)
  #### SummarizedExperiment ####
  } else if(class %in% c("summarizedexperiment","se")){
    obj <- to_se(obj = base_obj,
                 as_sce = FALSE,
                 verbose = verbose,
                 ...)
  } else if(class %in% c("singellcellexperiment","sce")){
    obj <- to_se(obj = base_obj,
                 as_sce = TRUE,
                 verbose = verbose,
                 ...)
  #### SingleCellExperiment ####
  } else if(class %in% c("seurat","s")){
      obj <- base_obj
  #### hdf5se ####
  } else if(class %in% c("hdf5se")){
    obj <- seurat_to_hdf5se(obj = base_obj,
                            verbose = verbose,
                            save_path = save_path,
                            ...)
    #### SeuratDisk ####
  } else if(class %in% c("h5seurat","h5s")){
    obj <- to_seurat(obj = obj,
                     as_h5seurat = TRUE,
                     save_path = paste0(save_path,".h5Seurat"),
                     verbose = verbose,
                     ...)
  #### CellDataSet ####
  } else if(class %in% c("celldataset","cds")){
    obj <- seurat_to_cds(obj = base_obj,
                         verbose = verbose)
  #### loom ####
  } else if(class %in% c("loom")){
    obj <- to_loom(obj = obj,
                   save_path = save_path,
                   verbose = verbose,
                   ...)
  #### LoomExperiment ####
  } else if(class %in% c("LoomExperiment","le")){
    obj <- seurat_to_le(obj =  base_obj,
                        as_scle = FALSE,
                        verbose = verbose)
  #### SingleCellLoomExperiment ####
  } else if(class %in% c("SingleCellLoomExperiment","scle")){
    obj <- seurat_to_le(obj =  base_obj,
                        as_scle = TRUE,
                        verbose = verbose)
  #### hdf5array ####
  } else if(class %in% c("hdf5array","h5")){
    obj <- seurat_to_hdf5se(obj = obj,
                            save_path = save_path,
                            verbose = verbose)
  #### anndata ####
  } else if(class %in% c("anndata","ad")){
      obj <- example_anndata(obj = base_obj,
                             save_path = save_path,
                             verbose = verbose,
                             ...)
  #### EWCE CelltypeDataset ####
  } else if(class %in% c("celltypedataset","celltypedata","ctd","ewce")){
    obj <- ewceData::ctd()
  #### list ####
  } else if(class %in% c("list","l")){
    obj <- example_list(obj = base_obj,
                        verbose = verbose,
                        ...)
  #### list paths ####
  }  else if(class %in% c("list_paths","lp")){
    obj <- example_list(obj = base_obj,
                        as_paths = TRUE,
                        verbose = verbose)
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
