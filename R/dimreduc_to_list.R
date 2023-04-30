dimreduc_to_list <- function(r){
  res <- list()
  slots <- methods::slotNames(r)
  if("cell.embeddings" %in% slots) {
    res[["obsm"]] <- r@cell.embeddings
  }
  if("feature.loadings" %in% slots) {
    res[["varm"]] <- r@feature.loadings
  }
  if("feature.loadings.projected" %in% slots) {
    res[["varm_projected"]] <- r@feature.loadings.projected
  }
  return(res)
}
