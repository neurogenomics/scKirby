#' Get unstructured metadata
#'
#' Extract unstructured metadata from any single-cell object.
#' @inheritParams converters
#' @inheritParams get_n_elements
#' @returns Named list of unstructured objects.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' uns <- get_uns(obj)
get_uns <- function(obj,
                    n=NULL,
                    verbose=TRUE){
  # devoptera::args2vars(get_uns)

  #### loom ####
  if(is_class(obj,"loom")){
    uns <- SeuratDisk::as.list(obj[["attrs"]])
  #### SummarizedExperiment ####
  } else if(is_class(obj,"se")){
    uns <- obj@metadata
  #### Seurat ####
  } else if(is_class(obj,"seurat")){
    ## Seurat V1
    if(methods::is(obj,"seurat")){
      uns <- obj@misc
      ## Seurat V2+
    } else {
      uns <- list(misc=obj@misc,
                  var_features=lapply(obj@assays,function(a){a@var.features})
                  )
    }
    #### h5Seurat ####
  } else if(is_class(obj,"h5seurat")){
    uns <- SeuratDisk::as.list(obj[["misc"]])
    #### anndata ####
  } else if(is_class(obj,"anndata")){
    uns <- obj$uns
    #### CellDataSet ####
  } else if(is_class(obj,"cds")){
    uns <- lapply(stats::setNames(methods::slotNames(obj@experimentData),
                                  methods::slotNames(obj@experimentData)),
                  function(x){
                    methods::layer(obj@experimentData,x)
                  })
    #### list ####
  } else if(is_class(obj,"list")){
    #### File path ####
    if(is.character(obj$uns)){
      uns <- read_data(path = obj$uns,
                       verbose = verbose,
                       as_sparse = FALSE)
    } else {
      uns <- obj$uns
    }
    #### OTHER ####
  } else {
    messager("Unable to get `uns` from object.",v=verbose)
  }
  #### Return as a named list (1 per assay), unless there's only 1 assay ####
  uns <- get_n_elements(l = uns,
                        n = n,
                        verbose = verbose)
  #### Return ####
  return(uns)
}
