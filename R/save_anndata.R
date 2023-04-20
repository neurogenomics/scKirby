save_anndata <- function(obj,
                         save_path,
                         method = c("anndata","zellkonverter"),
                         compression = "gzip",
                         compression_opts = NULL,
                         as_dense = list(),
                         conda_env="r-reticulate",
                         verbose=TRUE,
                         ...){
  obj <- example_obj("ad")
  method <- method[1]
  if(!is_class(obj,"anndata")){
    stopper("obj must be an anndata object.")
  }
  messager("+ Saving anndata:",save_path,v=verbose)
  if(method=="anndata"){
    echoconda::activate_env(conda_env = conda_env,
                            method = "reticulate",
                            verbose = verbose)
    anndata::write_h5ad(anndata = obj,
                        filename = save_path,
                        compression = compression,
                        compression_opts = compression_opts,
                        as_dense = as_dense)
  } else {
    ### Has to convert to SCE first ####
    obj <- zellkonverter::AnnData2SCE(obj,
                                      verbose = verbose)
    zellkonverter::writeH5AD(sce = obj,
                             file = save_path,
                             verbose = verbose,
                             ...)
  }

}
