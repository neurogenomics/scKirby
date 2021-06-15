
convert_to_sce <- function(object,
                           verbose=T){
  messager("Converting formats:",v=verbose)
  core_allocation <- assign_cores(worker_cores = .90)

  ewce_classes <- c("EWCE_list","SCE_list")
  matrix_classes <- c("data.table","data.frame","tbl_df","tbl","matrix","Matrix","array","DelayedArray","DelayedMatrix",names(getClass("Matrix")@subclasses)) # more than 40 ..
  loom_classes <-  c("loom","H5File","H5RefClass")
  sce_classes <- c("SingleCellLoomExperiment","SingleCellExperiment","SummarizedExperiment")
  anndata_classes <- c("AnnData","AnnDataR6","AnnDataR6R6")
  seurat_classes <- c("Seurat")
  supported_classes <- c(matrix_classes, loom_classes, sce_classes, anndata_classes, seurat_classes, ewce_classes)
  supported_classes_print <- c("matrix (any subclass)",loom_classes, sce_classes, anndata_classes, seurat_classes, ewce_classes)

  #### EWCE_list class ####
  if(class(object)[1]=="list" & all(c("exp","annot") %in% names(object)) ){
    class(object) <- "EWCE_list"
  }
  #### CTD class ####
  if(class(object)[1]=="list" & length(object)>1){
    if(all(c("mean_exp","specificity") %in% names(object[[1]]) )){
      messager("+ CTD ==> SingleCellExperiment",v=verbose)
      # object <- readRDS("~/Desktop/model_celltype_conservation/processed_data/EWCE/CTD_list.rds")[[1]]
      ctd <- object
      #### Name CTD levels ####
      if(is.null(names(ctd))){
        names(ctd) <- paste0("level_",1:length(ctd))
      } else {
        names(ctd) <- names(ctd)
      }
      sce_list <- lapply(names(ctd), function(lvl){
        message("Converting level: ",lvl)
        ctd_lvl <- ctd[[lvl]]
        #### Use matrices that are present ###
        matrix_list <- list()
        for(mtx_name in c("mean_exp","median_exp",
                          "specificity","median_specificity","specificity_quantiles")){
          if(mtx_name %in% names(ctd_lvl)){ matrix_list[[mtx_name]] <- DelayedArray::DelayedArray(as(as(ctd_lvl[[mtx_name]], "matrix"), "sparseMatrix")) }
        }
        sce <- SingleCellExperiment::SingleCellExperiment(
          assays      =  matrix_list,
          colData     =  data.frame(colnames(matrix_list[[1]])) %>% `colnames<-`(lvl),
          rowData     =  data.frame(gene=row.names(matrix_list[[1]]), row.names = row.names(matrix_list[[1]]))
        )
        sce <- check_sce_rownames(sce, verbose = verbose)
      }) %>% `names<-`(names(ctd))
      ## "SCE_list" class messes up other functions that expect class "list"
      # class(sce_list) <- "SCE_list"
      return(sce_list)
    }
  }

  #### Check if class is supported ####
  if(!class(object)[1] %in% supported_classes){
    stop("Unsupported class detected: ",class(object),"\n\n",
         "Data object must be at least one of the following classes:\n\n",
         paste(supported_classes_print, collapse = ", "))
  }
  ### EWCE_list ####
  if(class(object)[1] == "EWCE_list"){
    messager("+ EWCE_list ==> SingleCellExperiment",v=verbose)
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays      = list(raw = DelayedArray::DelayedArray(as(as.matrix(object$exp), "sparseMatrix"))),
      colData     =  object$annot,
      rowData     =  row.names(object$exp)
    )
    sce <- check_sce_rownames(sce, verbose = verbose)
    return(sce)
  }

  #### Matrices ####
  if(class(object)[1] %in% matrix_classes){
    messager("+ Matrix ==> SingleCellExperiment",v=verbose)
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays      = list(raw = DelayedArray::DelayedArray(as(as.matrix(object), "sparseMatrix"))),
    )
    sce <- check_sce_rownames(sce, verbose = verbose)
    return(sce)
  }
  #### Seurat ####
  if(class(object)[1] %in% seurat_classes){
    messager("+ Seurat ==> SingleCellExperiment",v=verbose)
    # object <- Seurat::pbmc_small ## example
    ### Seurat::as.SingleCellExperiment() is currently broken
    # sce <- Seurat::as.SingleCellExperiment(object)

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
  #### AnnData ####
  if(class(object)[1] %in% anndata_classes){
    messager("+ AnnData ==> SingleCellExperiment",v=verbose)
    # sceasy::convertFormat(obj = object, from = "anndata", to="seurat",
    #                       outFile = "~/Desktop/tmp.rds")
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays      = list(raw = DelayedArray::DelayedArray(as( Matrix::t(object$X), "sparseMatrix"))),
      colData     = object$obs,
      rowData     = object$var
    )
    sce <- check_sce_rownames(sce, verbose = verbose)
    return(sce)
  }
  if(class(object)[1] %in% loom_classes){
    messager("+ loom ==> SingleCellExperiment",v=verbose)
    # Import as a Seurat object first for convenience
    object <- SeuratDisk::LoadLoom(file = object$filename)
    assay1 <- names(object@assays)
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays      = list(raw = DelayedArray::DelayedArray(as(object@assays[[assay1]]@counts, "sparseMatrix"))),
      colData     = object@meta.data,
      rowData     = object@assays[[assay1]]@meta.features
    )
    sce <- check_sce_rownames(sce, verbose = verbose)
    return(sce)
  }
  #### SingleCellExperiment/SummarizedExperiment ####
  if(class(object)[1] %in% sce_classes){
    messager("+ == SummarizedExperiment",v=verbose)
    sce <- object
    sce <- check_sce_rownames(sce, verbose = verbose)
    return(sce)
  }
}
