#' Class dictionary
#'
#' Single-cell object class dictionary.
#' @param group Name of one or more single-cell group of
#'  single-cell object classes.
#' @returns Named list of object classes.
#'
#' @export
#' @examples
#' classes <- dict_class()
dict_class <- function(group=NULL){

  classes <- list(
    ewce = c("EWCElist","EWCE_list","SCElist","SCE_list",
             "CellTypeDataset"),
    matrix = c("data.table","data.frame","tbl_df","tbl",
               "matrix","Matrix","array","DelayedArray","DelayedMatrix",
               names(methods::getClass("Matrix")@subclasses),
               "Graph"
               ),
    matrix_strict=c("matrix","Matrix","array","DelayedArray","DelayedMatrix",
                    names(methods::getClass("Matrix")@subclasses)
                    ),
    sparse_matrix=names(methods::getClass("sparseMatrix")@subclasses),
    delayed_array=c("DelayedArray","DelayedMatrix"),
    loom = c("loom","H5File","H5RefClass"),
    se = c("SummarizedExperiment",
           "SingleCellLoomExperiment",
           "LoomExperiment",
           "SingleCellExperiment",
           "sce","se","le","scle"),
    hdf5se = c("HDF5SummarizedExperiment"),
    anndata = c("AnnData","AnnDataR6","AnnDataR6R6",
                "anndata._core.anndata.AnnData",
                "HDF5AnnData","InMemoryAnnData"),
    seurat = c("Seurat","SeuratObject","seurat"),
    h5seurat = c("h5Seurat","scdisk"),
    cds = c("ExpressionSet","CellDataSet","monocle","monocle3"),
    list = "list"
  )
  classes$supported <- unlist(classes, use.names = FALSE)
  classes$supported_print <- c(
    "matrix/data.frame (any subclass)",
    unlist(classes[names(classes)!="matrix"], use.names = FALSE)
  )
  if(is.null(group)){
    return(classes)
  } else{
    group_valid <- group[group %in% names(classes)]
    select_classes <- classes[group_valid]
    if(length(select_classes)==1){
      return(select_classes[[1]])
    } else {
      return(select_classes)
    }
  }
}
