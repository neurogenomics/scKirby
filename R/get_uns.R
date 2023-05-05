#' Get unstructured metadata
#'
#' Extract unstructuerd metadata from any single-cell object.
#' @export
#' @examples
#' obj <- example_obj("s")
#' uns <- get_uns(obj)
get_uns <- function(obj,
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
                    methods::slot(obj@experimentData,x)
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
  #### Return ####
  return(uns)
}
