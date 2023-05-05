#' Convert: \code{Seurat} ==> \code{AnnData}
#'
#' @param reimport Save and re-import the \code{AnnData} object into R to ensure
#' all data has been converted from Python-native to R-native objects
#' (e.g. pandas data.frames vs. R data.frames).
#' @inheritDotParams anndata::write_h5ad
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- seurat_to_anndata(obj)
seurat_to_anndata <- function(obj,
                              reimport = TRUE,
                              save_path = tempfile(fileext = ".h5ad"),
                              verbose = TRUE,
                              ...){

  messager("+ Seurat ==> AnnData",v=verbose)
  activate_conda(verbose=verbose)
  #### Convert ####
  adat <- sceasy::convertFormat(obj = obj,
                                from = "seurat",
                                to = "anndata",
                                outFile = NULL)
  ##### Save ####
  if(!is.null(save_path)){
    save_anndata(adat = adat,
                 save_path = save_path,
                 verbose = verbose)
  }
  #### Reimport ####
  if(isTRUE(reimport)){
    adat <- reimport_anndata(adat = adat,
                             save_path = save_path,
                             verbose = verbose)
  }
  #### Return ###
  return(adat)
}
