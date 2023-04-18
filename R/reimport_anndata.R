reimport_anndata <- function(obj,
                             save_path = tempfile(fileext = ".h5ad"),
                             verbose = TRUE,
                             ...){
  dir.create(dirname(save_path),showWarnings = FALSE, recursive = TRUE)
  # IMPORTANT! python cannot interpret "~"
  save_path <- path.expand(save_path)
  save_anndata(obj = obj,
               save_path = save_path,
               verbose = verbose,
               ...)
  obj <- anndata::read_h5ad(save_path)
  return(obj)
}
