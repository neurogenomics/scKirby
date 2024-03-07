#' Convert: \code{SummarizedExperiment} ==> \code{loom}
#'
#' @describeIn converters Convert
#' \link[SummarizedExperiment]{SummarizedExperiment} to
#' \link[SeuratDisk]{loom}
#' @inheritParams converters
#' @inheritDotParams SeuratDisk::SaveLoom
#' @returns description
#' @export
#' @examples
#' obj <- example_obj("se")
#' obj2 <- se_to_loom(obj)
se_to_loom <- function(obj,
                       save_path = file.path(tempdir(),
                                             "scKirby.loom"),
                       verbose=TRUE,
                       ...){
  messager_to()
  dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
  obj2 <- SeuratDisk::SaveLoom(object = obj,
                               filename = save_path,
                               verbose = verbose,
                               ...)
  return(obj2)
}
