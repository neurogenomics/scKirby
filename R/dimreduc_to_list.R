#' Convert: \code{AnnData} ==> \code{list}
#'
#' @inheritParams converters
#' @export
#' @examples
#' obj <- example_obj("anndata")
#' obj2 <- anndata_to_list(obj)
dimreduc_to_list <- function(obj){
  res <- list()
  slots <- methods::slotNames(obj)
  if("cell.embeddings" %in% slots) {
    res[["obsm"]] <- obj@cell.embeddings
  }
  if("feature.loadings" %in% slots) {
    res[["varm"]] <- obj@feature.loadings
  }
  if("feature.loadings.projected" %in% slots) {
    res[["varm_projected"]] <- obj@feature.loadings.projected
  }
  return(res)
}
