reimport_anndata <- function(adat,
                             save_path = tempfile(fileext = ".h5ad"),
                             verbose = TRUE){

  if(!is.null(save_path) &&
     file.exists(save_path)){
    messager("Reimporting adat.",v=verbose)
    #### Read back in ####
    adat <- read_anndata(path = save_path,
                         verbose = FALSE)
  } else {
    messager("Unable to reimport adat from save_path:",save_path,v=verbose)
  }
  #### Return ####
  return(adat)
}
