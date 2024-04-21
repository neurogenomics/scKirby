#' Get data
#'
#' Extract expression matrix metadata from any single-cell object.
#' @param as_sparse Convert to a \link[Matrix]{sparseMatrix}.
#' @inheritParams converters
#' @inheritParams get_n_elements
#' @inheritParams to_sparse
#' @inheritParams SeuratObject::CreateSeuratObject
#' @inheritParams SeuratObject::Assays
#' @returns A named list of matrices.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' X <- get_x(obj)
get_x <- function(obj,
                  transpose=FALSE,
                  n=NULL,
                  simplify=TRUE,
                  assay=NULL,
                  slot=NULL,
                  as_sparse=FALSE,
                  verbose=TRUE){
  # devoptera::args2vars(get_x)

  #### matrix ####
  if(is_class(obj,"matrix")){
    data <- to_sparse(obj = obj,
                      verbose = verbose)
    #### loom ####
  } else if(is_class(obj,"matrix_list")){
    data <- obj
  }else if(is_class(obj,"loom")){
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
    }
  #### SummarizedExperiment ####
  } else if(is_class(obj,"se")){
    data <- obj@assays@data
  #### Seurat ####
  } else if(is_class(obj,"seurat")){
    ## Seurat V1
    if(methods::is(obj,"seurat")){
      data <- list(RNA.counts=obj@raw.data,
                   RNA.data=obj@data,
                   RNA.scale.data=obj@scale.data)
    ## Seurat V2+
    } else {
      assays <- obj@assays
      if(!is.null(assay)) assays[assays %in% assay]
      data <- lapply(assays,function(a){
        slots <- c("counts","data","scale.data")
        slots <- slots[sapply(slots,function(s){methods::.hasSlot(a,s)})]
        if(!is.null(slot)) slots <- slots[slots %in% slot]
        lapply(stats::setNames(slots,slots), function(s){
          methods::slot(a,s)
        })
      }) |> unlist(recursive = FALSE)
    }
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
    #### Read file ####
    if(is.character(obj$data)){
      data <- read_data(path = obj$data,
                        verbose = verbose)
    } else {
      data <- obj$data
    }
    #### OTHER ####
  } else {
    stopper("Unable to get data from object.")
  }
  #### Conver to sparse matrices ####
  if(isTRUE(as_sparse)){
    data <- lapply(data, to_sparse, verbose=verbose)
  } else {
    data <- lapply(data, as.matrix)
  }
  #### Tranpose each matrix ####
  if(isTRUE(transpose)){
    messager("Transposing data.",v=verbose)
    data <- lapply(data, Matrix::t)
  }
  #### Return as a named list (1 per assay), unless there's only 1 assay ####
  data <- get_n_elements(l = data,
                         n = n,
                         simplify = simplify,
                         verbose = verbose)
  #### Return ####
  return(data)
}
