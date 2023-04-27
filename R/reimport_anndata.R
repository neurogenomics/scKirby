reimport_anndata <- function(adat,
                             save_path = tempfile(fileext = ".h5ad"),
                             verbose = TRUE){

  if(!is.null(save_path) &&
     file.exists(save_path)){
    messager("Reimporting adat.",v=verbose)
    #### Read back in ####
    adat <- anndata::read_h5ad(save_path)
    #### Add filename to adat ####
    adat$filename <- save_path
  } else {
    messager("Unable to reimport adat from save_path:",save_path,v=verbose)
  }
  #### Return ####
  return(adat)
}
