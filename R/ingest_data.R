#' Import and standardize scRNAseq data across different formats
#'
#' Automatically infers data format of scRNAseq object,
#' or a path to that object.
#' It then uses the appropriate functions to import that data and convert
#'  it to a \pkg{SingleCellExperiment},
#'  which is recognized by other \pkg{EWCE} functions.
#'
#' \href{https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html}{SeuratDisk}
#' \href{https://github.com/rcannood/anndata}{anndata (R)}
#' \href{https://anndata.readthedocs.io/en/latest/}{anndata (python)}
#' \href{https://satijalab.org/loomR/loomR_tutorial.html}{loomR}
#' \href{https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html}{SingleCellExperiment}
#' \href{https://petehaitch.github.io/BioC2020_DelayedArray_workshop/articles/Effectively_using_the_DelayedArray_framework_for_users.html}{DelayedArray workshop}
#' \href{https://theislab.github.io/zellkonverter/articles/zellkonverter.html}{zellkonverter}
#' @param obj A single-cell data object, or a path to saved single-cell data.
#' @param input_type Format of \code{obj}. By default, the type will be inferred.
#' @param output_class Format to convert \code{obj} to.
#' @param custom_reader Custom function to read \code{obj} into R.
#' @param save_path Path to save the converted \code{obj}.
#' @param overwrite If a file of the same name exists, overwrite it.
#' @param return_save_path If \code{TRUE}, a list with both the converted object
#'  and the saved file path will be returned
#'  (instead of just the converted object).
#' @inheritParams converters
#' @inheritDotParams read_data
#' @returns Converted single-cell object.
#'
#' @export
#' @examples
#' \dontrun{
#' #### Input: Seurat object ####
#' seurat <- example_obj("seurat")
#' sce1 <- ingest_data(obj=seurat)
#'
#' #### Input: Matrix object ####
#' sce2 <- ingest_data(obj=seurat@assays$RNA@counts)
#'
#' #### Input hdf5se object ####
#' hdf5se <- example_obj("hdf5se")
#' sce3 <- ingest_data(obj=hdf5se)
#'
#' #### Input anndata object ####
#' ad <- example_obj("anndata")
#' sce4 <- ingest_data(obj=ad)
#'
#' #### Input anndata path ####
#' sce5 <- ingest_data(obj = ad$filename)
#'
#' #### Input: h5Seurat object ####
#' h5seurat <- example_obj("h5seurat")
#' sce6 <- ingest_data(obj=h5seurat)
#'
#' #### Input: h5Seurat path ####
#' sce7 <- ingest_data(obj=h5seurat$path)
#'
#' #### Input: loom object ####
#' library(Seurat)
#' loom <- example_obj("loom")
#' sce8 <- ingest_data(obj=loom)
#'
#' #### Input: EWCE object ####
#' \dontrun{
#' ewce <- ewceData::cortex_mrna()
#' sce10 <- ingest_data(obj=ewce)
#' }
#' }
ingest_data <- function(obj,
                        input_type = "guess",
                        output_class = c("SingleCellExperiment",
                                         "Seurat",
                                         "CellDataSet",
                                         "list"),
                        custom_reader = NULL,
                        save_path = NULL,
                        overwrite = FALSE,
                        return_save_path = FALSE,
                        verbose = TRUE,
                        ...){
  # devoptera::args2vars(ingest_data)

  #### Select output type ####
  output_class <- dict_output(output_class = output_class)
  #### Read ####
  obj <- read_data(path=obj,
                   filetype=input_type,
                   custom_reader=custom_reader,
                   verbose=verbose,
                   ...)
  #### Convert ####
  cdict <- dict_class()
  #### to SummarizedExperiment / SingleCellExperiment ####
  if(output_class %in% tolower(cdict$se)){
      obj_out <- to_se(obj = obj,
                       as_sce = output_class=="singlecellexperiment",
                       verbose = verbose)
  #### to Seurat ####
  } else if(output_class %in% tolower(c(cdict$seurat,cdict$h5seurat))){
      obj_out <- to_seurat(obj = obj,
                           as_h5seurat = is_filetype(output_class,"h5seurat"),
                           verbose = verbose)
  #### to Loom ####
  } else if(output_class %in% tolower(cdict$loom)){
    obj_out <- to_loom(obj = obj,
                       save_path = save_path,
                       verbose = verbose)
  #### to AnnData ####
  }else if(output_class %in% tolower(cdict$anndata)){
    obj_out <- to_anndata(obj = obj,
                          save_path = save_path,
                          verbose = verbose)
  }
  #### Save ####
  if(!is.null(save_path)){
    save_path <- save_data(obj=obj_out,
                           filetype=output_class,
                           save_path=save_path,
                           overwrite=overwrite,
                           verbose=verbose)
  }
  #### Return ####
  if(isTRUE(return_save_path)) {
    return(list(save_path=save_path,
                obj=obj_out))
  } else {
    return(obj_out)
  }
}



