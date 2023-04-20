save_anndata <- function(adat,
                         save_path,
                         method = c("anndata","zellkonverter"),
                         compression = "gzip",
                         compression_opts = NULL,
                         as_dense = list(),
                         conda_env="r-reticulate",
                         verbose=TRUE,
                         ...){
  # devoptera::args2vars(save_anndata)

  method <- method[1]
  if(!is_class(adat,"anndata")){
    stopper("obj must be an anndata object.")
  }
  # IMPORTANT! python cannot interpret "~"
  save_path <- path.expand(save_path)
  #### Make sure directory exists ####
  dir.create(dirname(save_path),showWarnings = FALSE, recursive = TRUE)
  #### Save anndata ####
  messager("+ Saving anndata -->",save_path,v=verbose)
  if(method=="anndata"){
    echoconda::activate_env(conda_env = conda_env,
                            method = "reticulate",
                            verbose = verbose)
    adat$write_h5ad(filename = save_path,
                    compression = compression,
                    compression_opts = compression_opts,
                    as_dense = as_dense)
  } else {
    ### Has to convert to SCE first ####
    obj <- zellkonverter::AnnData2SCE(adat,
                                      verbose = verbose)
    zellkonverter::writeH5AD(sce = obj,
                             file = save_path,
                             verbose = verbose,
                             ...)
  }
  #### Add filename to adat ####
  if(!is.null(save_path) &&
     file.exists(save_path)){
    adat$filename <- save_path
  }
}
