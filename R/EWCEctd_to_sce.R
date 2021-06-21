
#' Convert: \code{EWCEctd} ==> \code{SingleCellExperiment}
#'
#' @examples
#' library(scKirby)
#' data("example_ctd")
#' sce <- EWCEctd_to_sce(example_ctd)
EWCEctd_to_sce <- function(object,
                           verbose=T){
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
