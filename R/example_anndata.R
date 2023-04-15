#' Example \code{anndata}
#'
#' \href{https://github.com/dynverse/anndata/issues/1}{
#' anndata R package has a major bug currently.}
#' Thus, \code{add_filename} must be kept \code{False} until this is fixed.
#' @source
#' \code{
#' #### Setup conda env ####
#' ## Option 1: Point to where anndata is installed (or should be installed)
#' conda_dir <- dirname(dirname(reticulate::conda_list()[1,]$python))
#' reticulate::use_condaenv(condaenv = conda_dir)
#' reticulate::conda_install(conda = conda_dir, packages = "loompy", pip = T)
#' ## Option 2: Just run anndata::install_anndata() and will install via miniconda
#' anndata::install_anndata(method = "conda", conda=conda_dir)
#' #### Usage ####
#' adata <- example_anndata()
#' }
#' @inheritParams echoconda::activate_env
#'
#' @keywords internal
#' @importFrom echoconda activate_env
example_anndata <- function(obj = SeuratObject::pbmc_small,
                            return_path=FALSE,
                            backed=c("r","r+"),
                            overwrite=TRUE,
                            add_filename=FALSE,
                            save_path = file.path(tempdir,"example.h5ad"),
                            conda_env = "r-reticulate",
                            verbose=TRUE){

  echoconda::activate_env(conda_env = conda_env,
                          method = "reticulate",
                          verbose = verbose)

  backed <- backed[1]
  if(file.exists(save_path) &&
     isFALSE(overwrite)){
    messager("+ Importing existing anndata object", save_path, v=verbose)
    adat <- anndata::read_h5ad(filename = save_path,
                               backed = backed)
  } else {
    messager("+ Creating new anndata object:",save_path, v=verbose)
    t_mat <-  Matrix::t(Seurat::GetAssayData(obj))
    adat <- anndata::AnnData(
      X = t_mat,
      obs = obj@meta.data[rownames(t_mat),],
      var = Seurat::GetAssay(obj)@meta.features[colnames(t_mat),]
      )
    if(is.null(save_path)) return (adat)
    if(isTRUE(overwrite)) out <- suppressWarnings(file.remove(save_path))
    out <- adat$write_h5ad(filename = save_path)
    if(isTRUE(add_filename)) adat$filename <- save_path
  }
  if(isTRUE(return_path)) return(save_path) else return(adat)
}
