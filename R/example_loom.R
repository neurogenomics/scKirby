

## Certain file types (loom, anndata) can't be stored a .rda
## because they need to be set up for on-disk reading.


#' Example \code{loom}
#'
#' @examples
#' \dontrun{
#' #### Usage ####
#' library(scKirby)
#' loom <- example_loom()
#' @export
example_loom <- function(save_dir=tempdir()){
  library(Seurat) # MUST be loaded
  library(SeuratDisk)
  loom <- SeuratDisk::as.loom(scKirby::example_seurat,
                              overwrite=T,
                              filename=file.path(save_dir,"pbmc_small.loom"))
  return(loom)
}



