

class_dictionary <- function(){
  ewce_classes <- c("EWCElist","EWCE_list","SCElist","SCE_list")
  matrix_classes <- c("data.table","data.frame","tbl_df","tbl","matrix","Matrix","array","DelayedArray","DelayedMatrix",names(getClass("Matrix")@subclasses)) # more than 40 ..
  loom_classes <-  c("loom","H5File","H5RefClass")
  sce_classes <- c("SingleCellLoomExperiment","SingleCellExperiment","SummarizedExperiment","sce")
  hdf5se_classes <- c("HDF5SummarizedExperiment")
  anndata_classes <- c("AnnData","AnnDataR6","AnnDataR6R6")
  seurat_classes <- c("Seurat","SeuratObject")
  seurath5_classes <- c("h5Seurat","scdisk")
  cds_classes <- c("ExpressionSet","CellDataSet","monocle","monocle3")
  #### Aggregate ####
  supported_classes <- c(matrix_classes, loom_classes, sce_classes, hdf5se_classes, anndata_classes, seurat_classes, seurath5_classes, ewce_classes, cds_classes)
  supported_classes_print <- c("matrix (any subclass)",loom_classes, sce_classes, hdf5se_classes, anndata_classes, seurat_classes, ewce_classes, cds_classes)

  return(list(ewce=ewce_classes,
              matrix=matrix_classes,
              loom=loom_classes,
              sce=sce_classes,
              hdf5se=hdf5se_classes,
              anndata=anndata_classes,
              seurat=seurat_classes,
              h5seurat=seurath5_classes,
              cds=cds_classes,

              supported_classes=supported_classes,
              supported_classes_print=supported_classes_print
  ))
}


