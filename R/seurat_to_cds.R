#' Convert: \code{Seurat} ==> \code{CellDataSet}
#'
#' NOTE: \pkg{Seurat} does not take into account the current \pkg{monocle3}
#' format: \link[monocle3]{cell_data_set}. Instead, it only converts to the old
#' \link[monocle]{CellDataSet} format.
#' @inheritParams converters
#' @returns A \link[monocle3]{cell_data_set} or
#' \link[monocle]{CellDataSet} object.
#'
#' @export
#' @importFrom Seurat as.CellDataSet
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- seurat_to_cds(obj)
seurat_to_cds <- function(obj,
                          version=c("monocle3","monocle"),
                          verbose=TRUE,
                          ...){

  version <- tolower(version)[1]
  #### CellDataSet ####
  if(version=="monocle"){
    Seurat::as.CellDataSet(obj)
  #### cell_data_set ####
  } else if(version=="monocle3"){
    l <- to_list(obj = obj, verbose = verbose)
    list_to_cds(obj = l,
                verbose = verbose,
                ...)
  } else {
    stopper("verision must be one of",
            paste("\n -",shQuote(eval(formals(seurat_to_cds)$version)),
                  collapse = "")
            )
  }
}
