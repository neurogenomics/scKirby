#' Get data
#'
#' Extract expression matrix metadata from any single-cell object.
#' @export
#' @examples
#' obj <- example_obj("loom")
#' X <- get_data(obj)
get_data <- function(obj,
                     transpose=TRUE,
                     simplify=FALSE,
                     verbose=TRUE){
  # devoptera::args2vars(get_data)

  #### matrix ####
  if(is_class(obj,"matrix")){
    data <- to_sparse(obj = obj,
                      verbose = verbose)
    #### loom ####
  } else if(is_class(obj,"loom")){
    data <- list()
    if("matrix" %in% names(obj)){
      data$X <- Seurat::as.sparse(obj[["matrix"]])
      rownames(data$X) <- get_obs(obj)[[1]]
      colnames(data$X) <- get_var(obj)[[1]]
    }
    if("layers" %in% names(obj)){
      nms <- names(obj[["layers"]])
      for(nm in nms){
        data[[nm]] <- Seurat::as.sparse(obj[["layers"]][[nm]])
        rownames(data[[nm]]) <- get_obs(obj)[[1]]
        colnames(data[[nm]]) <- get_var(obj)[[1]]
      }
      if(isTRUE(transpose)){
        messager("Transposing data.",v=verbose)
        data <- lapply(data, Matrix::t)
      }
    }
  #### SummarizedExperiment ####
  } else if(is_class(obj,"se")){
    data <- obj@assays@data
  #### Seurat ####
  } else if(is_class(obj,"seurat")){
    data <- lapply(obj@assays,function(a){
      slots <- c("counts","data","scale.data")
      slots <- slots[sapply(slots,function(s){methods::.hasSlot(a,s)})]
      lapply(stats::setNames(slots,slots), function(s){
        methods::slot(a,s)
      })
    }) |> unlist(recursive = FALSE)
  #### h5Seurat ####
  } else if(is_class(obj,"h5seurat")){
    assays <- obj[["assays"]]$ls()$name
    data <- lapply(stats::setNames(assays,
                                   assays),
                  function(nm){
                    a <- obj[["assays"]][[nm]]
                    slots <- c("counts","data")
                    slots <- slots[slots %in% names(a)]
                    lapply(stats::setNames(slots,slots), function(s){
                      X <- SeuratObject::as.sparse(a[[s]])
                      colnames(X) <- rownames(get_obs(obj))
                      rownames(X) <- rownames(get_var(obj))
                      return(X)
                    })
                  }) |> unlist(recursive = FALSE)
  #### anndata ####
  } else if(is_class(obj,"anndata")){
    data <- obj$X
    #### CellDataSet ####
  } else if(is_class(obj,"cds")){
    data <- as.list(obj@assayData)
    #### list ####
  } else if(is_class(obj,"list")){
    data <- obj$data
    #### OTHER ####
  } else {
    stopper("Unable to get data from object.")
  }
  #### Return as a named list (1 per assay), unless there's only 1 assay ####
  if(isTRUE(simplify)){
    data <- simplify_list(l = data)
  }
  #### Return ####
  return(data)
}
