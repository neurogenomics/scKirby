
#' Convert: \code{EWCEctd} ==> \code{Seurat}
#'
#' @examples
#' library(scKirby)
#' data("example_ctd")
#' seurat <- EWCEctd_to_seurat(example_ctd)
EWCEctd_to_seurat <- function(object,
                              verbose=T){
  messager("+ Seurat ==> SingleCellExperiment",v=verbose)
  ctd <- object
  #### Name CTD levels ####
  if(is.null(names(ctd))){
    names(ctd) <- paste0("level_",1:length(ctd))
  } else {
    names(ctd) <- names(ctd)
  }
  seurat_list <- lapply(names(ctd), function(lvl){
    message("Converting level: ",lvl)
    ctd_lvl <- ctd[[lvl]]
    #### Use matrices that are present ###
    matrix_list <- list()
    for(mtx_name in c("mean_exp","median_exp",
                      "specificity","median_specificity","specificity_quantiles")){
      if(mtx_name %in% names(ctd_lvl)){ matrix_list[[mtx_name]] <- as(as(ctd_lvl[[mtx_name]], "matrix"), "sparseMatrix") }
    }
    seurat <- CreateSeuratObject(counts = matrix_list[[1]],
                                 meta.data = data.frame(colnames(matrix_list[[1]]),
                                                        row.names =colnames(matrix_list[[1]]) ) %>% `colnames<-`(lvl)
                                 )
    return(seurat)
  }) %>% `names<-`(names(ctd))
  return(seurat_list)
}
