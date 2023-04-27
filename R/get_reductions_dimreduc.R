get_reductions_dimreduc <- function(r){
  res <- list()
  slots <- methods::slotNames(r)
  if("cell.embeddings" %in% slots) {
    res[["embeddings"]] <- r@cell.embeddings
  }
  if("feature.loadings" %in% slots) {
    res[["loadings"]] <- r@feature.loadings
  }
  if("feature.loadings.projected" %in% slots) {
    res[["loadings_projected"]] <- r@feature.loadings.projected
  }
  return(res)
}
