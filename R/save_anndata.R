save_anndata <- function(adat,
                         save_path,
                         method = c("anndata","zellkonverter"),
                         compression = "gzip",
                         compression_opts = NULL,
                         as_dense = list(),
                         verbose=TRUE,
                         ...){
  # devoptera::args2vars(save_anndata)

  method <- method[1]
  if(!is_class(adat,"anndata")){
    stopper("obj must be an AnnData object.")
  }
  # IMPORTANT! python cannot interpret "~"
  save_path <- path.expand(save_path)
  #### Make sure directory exists ####
  dir.create(dirname(save_path),showWarnings = FALSE, recursive = TRUE)
  #### Save anndata ####
  messager("+ Saving AnnData -->",save_path,v=verbose)
  if(method=="anndata"){
    activate_conda(verbose = verbose)
    adat$write_h5ad(filename = save_path,
                    compression = compression,
                    compression_opts = compression_opts,
                    as_dense = as_dense)
  } else {
    #### Has to convert to SCE first ####
    obj <- zellkonverter::AnnData2SCE(adat,
                                      verbose = verbose)
    zellkonverter::writeH5AD(sce = obj,
                             file = save_path,
                             verbose = verbose,
                             ...)
  }
}
