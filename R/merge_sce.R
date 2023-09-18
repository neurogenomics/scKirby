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
#' @importFrom SingleCellExperiment colData
#' @importFrom EWCE merge_sce
#' @examples
#' obj <- example_obj("sce")
#' batch_names <- unique(obj$groups)
#' obj_list <- list(obj[,obj$groups=='g1'],
#'                  obj[,obj$groups=='g2'])
#' obj2 <- merge_sce(obj_list, batch_names = batch_names)
merge_sce <- function(obj_list,
                      method = "union",
                      cut_off_batch = 0.01,
                      cut_off_overall = 0.01,
                      use_assays = NULL,
                      colData_names = intersect,
                      batch_names = names(obj_list),
                      verbose = TRUE) {
  # devoptera::args2vars(merge_sce)

  if(is.function(colData_names)){
    colData_names <- lapply(obj_list, function(sce){
      colnames(SummarizedExperiment::colData(sce))
    }) |> do.call(what=colData_names)
  }
  EWCE::merge_sce(obj_list = obj_list,
                  method = method,
                  cut_off_batch = cut_off_batch,
                  cut_off_overall = cut_off_overall,
                  use_assays = use_assays,
                  colData_names = colData_names,
                  batch_names = batch_names,
                  verbose = verbose)
}
