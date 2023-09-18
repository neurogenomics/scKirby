#' Convert: to \code{loom}
#'
#' @describeIn converters
#' Convert any object to \link[loomR]{loom} format.
#' @inheritParams converters
#' @returns A \link[loomR]{loom} object.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- to_loom(obj)
to_loom <- function(obj,
                    save_path = file.path(tempdir(),
                                          "scKirby.loom"),
                    verbose=TRUE,
                    ...){
  if(is_class(obj,"loom")){
    return(obj)
  #### Seurat ####
  } else if(is_class(obj,"seurat")){
    if(!is.null(save_path)){
      obj <- SeuratDisk::as.loom(x = obj,
                                 filename = save_path,
                                 ...)
    } else {
      obj <- SeuratDisk::as.loom(x = obj,
                                 ...)
    }
  #### OTHER ####
  } else {
  obj <- ingest_data(obj = obj,
                     output_class = 'sce',
                     save_path = NULL,
                     verbose = verbose)
  obj2 <- se_to_loom(obj = obj,
                     save_path = save_path,
                     verbose = verbose)
  }
  return(obj2)
}
