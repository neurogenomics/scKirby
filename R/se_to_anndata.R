#' Convert: \code{SummarizedExperiment} ==> \code{AnnData}
#'
#' @param method Which R package to use for the conversion.
#' @param ... Additional arguments passed to the respective converter function.
#' @inheritParams converters
#' @inheritParams seurat_to_anndata
#' @returns \link[anndata]{AnnData} object.
#'
#' @export
#' @examples
#' obj <- example_obj("se")
#' obj2 <- se_to_anndata(obj)
se_to_anndata <- function(obj,
                          method = c("zellkonverter","sceasy"),
                          reimport = TRUE,
                          save_path = tempfile(fileext = ".h5ad"),
                          verbose=TRUE,
                          ...){

  method <- tolower(method[1])
  if(is_class(obj,"anndata")){
    messager("obj already in anndata format. Returning directly.",v=verbose)
    return(obj)
  }
  #### Ensure SCE format first ####
  obj <- se_to_sce(obj = obj,
                   verbose = verbose)
  #### Convert from HDF5Matrix ####
  ## Neither method can convert from this data type ###
  if(methods::is(obj@assays@data$counts,"HDF5Matrix")){
    obj@assays@data$counts <- as(obj@assays@data$counts,'sparseMatrix')
  }
  #### Convert ####
  messager("+ SummarizedExperiment ==> AnnData",v=verbose)
  if(method=="zellkonverter"){
    requireNamespace("zellkonverter")
    adat <- zellkonverter::SCE2AnnData(obj,
                               verbose = verbose,
                               ...)
  } else {
    adat <- sceasy::convertFormat(obj,
                          from = "sce",
                          to="anndata")
  }
  #### Reimport ####
  if(isTRUE(reimport)){
    adat <- reimport_anndata(adat = adat,
                             save_path = save_path,
                             ...)
  }
  return(adat)
}


