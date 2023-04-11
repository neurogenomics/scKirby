save_anndata <- function(obj,
                         save_path,
                         verbose=TRUE){
  if(!is_class(obj,"anndata")){
    stopper("obj must be an anndata object.")
  }
  messager("+ Saving anndata:",save_path,v=verbose)
  anndata::write_h5ad(anndata = obj,
                      filename = save_path)
}
