reimport_anndata <- function(adat,
                             save_path = tempfile(fileext = ".h5ad"),
                             verbose = TRUE,
                             ...){
  #### Read back in ####
  adat <- anndata::read_h5ad(save_path)
  #### Add filename to adat ####
  if(!is.null(save_path) &&
     file.exists(save_path)){
    adat$filename <- save_path
  }
  #### Return ####
  return(adat)
}
