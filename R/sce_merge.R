#' Merge a list of \code{SingleCellExperiment} objects
#'
#' Combine several \code{SingleCellExperiment} objects from
#' different batches/experiments.
#' Extracted from the
#'  \href{https://bioconductor.org/packages/release/bioc/html/scMerge.html}{
#'  scMerge package}.
#' @inheritParams EWCE::merge_sce
#' @returns A \code{SingleCellExperiment} object
#' with the list of SCE objects combined.
#'
#' @export
#' @import SingleCellExperiment
#' @examples
#' obj <- example_obj("sce")
#' batch_names<-unique(obj$batch)
#' sce_list<-list(obj[,obj$batch=='batch2'],
#'                obj[,obj$batch=='batch3'])
#' obj2 <- sce_merge(sce_list, batch_names = batch_names)
sce_merge <- function(sce_list,
                      method = "union",
                      cut_off_batch = 0.01,
                      cut_off_overall = 0.01,
                      use_assays = NULL,
                      colData_names = NULL,
                      batch_names = NULL,
                      verbose = TRUE) {
  EWCE::merge_sce(sce_list = sce_list,
                  method = method,
                  cut_off_batch = cut_off_batch,
                  cut_off_overall = cut_off_overall,
                  use_assays = use_assays,
                  colData_names = colData_names,
                  batch_names = batch_names,
                  verbose = verbose)
}
