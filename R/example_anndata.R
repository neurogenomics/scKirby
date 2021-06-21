
#' Example \code{anndata}
#'
#'
#' \href{https://github.com/dynverse/anndata/issues/1}{anndata R package has a major bug currently.}
#' Thus, \code{add_filename} must be kept \code{False} until this is fixed.
#'
#' @examples
#' \dontrun{
#' #### Setup conda env ####
#' ## Option 1: Point to where anndata is installed (or should be installed)
#' conda_dir <- dirname(dirname(reticulate::conda_list()[1,]$python))
#' reticulate::use_condaenv(condaenv = conda_dir)
#' reticulate::conda_install(conda = conda_dir, packages = "loompy", pip = T)
#' ## Option 2: Just run anndata::install_anndata() and will install via miniconda
#' anndata::install_anndata(method = "conda", conda=conda_dir)
#'
#' #### Usage ####
#' adata <- example_anndata()
#' }
#' @export
example_anndata <- function(return_path=F,
                            save_dir=tempdir(),
                            backed=c("r","r+"),
                            overwrite=T,
                            add_filename=F,
                            verbose=T){
  filename <- file.path(save_dir,"pbmc_small.h5ad")
  if(file.exists(filename) & overwrite==F){
    printer("+ Importing existing anndata object", filename, v=verbose)
    adata <- anndata::read_h5ad(filename = filename, backed = backed[1])
  }else {
    printer("+ Creating new anndata object:",filename, v=verbose)
    t_mat <-  SparseM::t(Seurat::GetAssayData(scKirby::example_seurat))
    adata <- anndata::AnnData(X = t_mat,
                              obs = scKirby::example_seurat@meta.data[rownames(t_mat),],
                              var = Seurat::GetAssay(scKirby::example_seurat)@meta.features[colnames(t_mat),]
                              )
    if(overwrite) out <- suppressWarnings(file.remove(filename))
    out <- adata$write_h5ad(filename = filename)
    if(add_filename) adata$filename <- filename
  }
  if(return_path) return(filename) else return(adata)
}
