#' Convert: to \code{CellDataSet}
#'
#' Convert any object to \code{CellDataSet} format.
#' @inheritParams converters
#' @returns A \code{CellDataSet} object.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- to_cds(obj)
to_cds <- function(obj,
                   save_path = file.path(tempdir(),
                                          "scKirby.loom"),
                   verbose=TRUE,
                  ...){

  if(is_class(obj,"cds")){
    return(obj)
    #### Seurat ####
  } else if(is_class(obj,"seurat")){
    obj2 <- seurat_to_cds(obj = obj,
                         verbose = verbose)
    #### OTHER ####
  } else {
    obj <- to_list(obj = obj,
                   verbose = verbose)
    obj2 <- list_to_cds(obj = obj,
                        verbose = verbose,
                        ...)
  }
  return(obj2)
}
