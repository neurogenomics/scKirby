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
#' @inheritParams anndata::read_h5ad
#'
#' @keywords internal
#' @importFrom echoconda activate_env
example_anndata <- function(obj = SeuratObject::pbmc_small,
                            save_path = file.path(tempdir(),"example.h5ad"),
                            backed = "r",
                            verbose = TRUE){
  # devoptera::args2vars(example_anndata)

  #### Import existing anndata file ####
  if(!is.null(save_path) &&
     file.exists(save_path)){
    adat <- read_anndata(path = save_path,
                         backed = backed,
                         verbose = verbose)
  } else {
    #### Create new  anndata file ####
    adat <- seurat_to_anndata(obj = obj,
                              save_path = save_path,
                              reimport = TRUE,
                              verbose = verbose)
    adat$filename <- save_path
  }
  return(adat)
}
