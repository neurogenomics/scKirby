#' Convert: \code{SummarizedExperiment} ==> \code{loom}
#'
#' @inheritParams converters
#' @export
#' @examples
#' obj <- example_obj("se")
#' obj2 <- se_to_loom(obj)
se_to_loom <- function(obj,
                       save_path = file.path(tempdir(),
                                             "scKirby.loom"),
                       verbose=TRUE){
  messager("+ SummarizedExperiment ==> loom",v=verbose)
  dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
  l <- to_list(obj)
  SeuratDisk::SaveLoom(object = obj)
  obj2 <- loomR::create(data = l$data,
                        filename = save_path)
  return(obj2)
}
