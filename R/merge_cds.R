#' Merge a list of \code{CellDataSet} objects
#'
#' Combine several \code{CellDataSet} objects from
#' different batches/experiments.
#' @returns A \code{CellDataSet} object.
#'
#' @export
#' @importFrom SummarizedExperiment assay assayNames
#' @examples
#' \dontrun{
#' obj <- example_obj("cds")
#' batch_names <- unique(obj$groups)
#' obj_list <- list(obj[,obj$groups=='g1'],
#'                  obj[,obj$groups=='g2'])
#' obj2 <- merge_cds(obj_list, batch_names = batch_names)
#' }
merge_cds <- function(obj_list,
                      by = c("assays","obs","var"),
                      verbose = TRUE) {
  # devoptera::args2vars(merge_cds)
  by <- tolower(by)[1]
  if(length(obj_list)==1){
    obj2 <- obj_list[[1]]
  } else {
    #### Merge by assays ####
    if(isTRUE(by=="assays")){
      messager("Merging CellDataSet objects as multiple assays.",v=verbose)
      obj2 <- obj_list[[1]]
      for(i in seq(2,length(obj_list))){
        SummarizedExperiment::assay(obj2, i = i) <- obj_list[[i]]
        SummarizedExperiment::assayNames(obj2)[i] <- names(obj_list)[i]
      }
    } else {
      #### Merge by obs ####
      messager("Merging CellDataSet objects as multiple obs (samples).",
               v=verbose)

      if(!requireNamespace("monocle3", quietly = TRUE)){
        stopper("monocle3 required for merging by obs.\n",
                "Install with: remotes::install_github('cole-trapnell-lab/monocle3')")
      }
      combine_cds <- utils::getFromNamespace("combine_cds", "monocle3")
      combine_cds(cds_list = obj_list)
    }

  }
  return(obj2)
}
