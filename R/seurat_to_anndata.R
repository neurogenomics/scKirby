#' Convert: \code{Seurat} ==> \code{AnnData}
#'
#' @param method R package to use when converting object.
#' @inheritParams converters
#' @inheritParams to_anndata
#' @inheritDotParams anndata::write_h5ad
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- seurat_to_anndata(obj, method="seuratdisk")
seurat_to_anndata <- function(obj,
                              reimport = TRUE,
                              save_path = tempfile(fileext = ".h5ad"),
                              method = c("sceasy","anndatar","seuratdisk"),
                              verbose = TRUE,
                              ...){

  messager_to()
  method <- tolower(method)[1]
  #### Convert ####
  if(!is.null(save_path)) save_path <- path.expand(save_path)
  # Method: anndatar
  if(method=="anndatar"){
    adat <- anndataR::from_Seurat(
      seurat_obj = obj,
      output_class = if(!is.null(save_path)){
        "HDF5AnnData"
      } else {
        "InMemoryAnnData"
      },
      file = save_path
    )
  } else if(method=="sceasy"){
    # Method: sceasy
    activate_conda(verbose=verbose)
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
  } else{
    messager("Saving ==>",save_path)
    h5Seurat_path <- gsub("\\.h5ad$",".h5Seurat",save_path)
    SeuratDisk::SaveH5Seurat(obj, filename = h5Seurat_path)
    SeuratDisk::Convert(h5Seurat_path, dest = "h5ad")
    adat <- anndata::read_h5ad(save_path)
  }
  #### Return ###
  return(adat)
}
