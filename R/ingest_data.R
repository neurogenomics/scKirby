#' Import and standardize scRNAseq data across different formats
#'
#' Automatically infers data format of scRNAseq object, or a path to that object.
#' It then uses the appropriate functions to import that data and convert it to a
#' \pkg{SingleCellExperiment}, which is recognized by other \pkg{EWCE} functions.
#'
#' @param obj A single-cell data object, or a path to saved single-cell data.
#' @param input_type Format of \code{obj}. By default, the type will be inferred.
#' @param output_class Format to convert \code{obj} to.
#' @param custom_reader Custom function to read \code{obj} into R.
#' @param save_dir Directory to save the converted \code{obj}.
#' @param filename Name to save the converted \code{obj}.
#' @param save_output Whether or not to save the converted \code{obj}.
#' @param overwrite If a file of the same name exists, overwrite it.
#' @param return_save_path If \code{TRUE}, a list with both the converted object
#'  and the saved file path will be returned (instead of just the converted object).
#' @param ... Additional arguments to be passed to \link[scKirby]{read_data}.
#' @source
#' \href{https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html}{SeuratDisk}
#' \href{https://github.com/rcannood/anndata}{anndata (R)}
#' \href{https://anndata.readthedocs.io/en/latest/}{anndata (python)}
#' \href{https://satijalab.org/loomR/loomR_tutorial.html}{loomR}
#' \href{https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html}{SingleCellExperiment}
#' \href{https://petehaitch.github.io/BioC2020_DelayedArray_workshop/articles/Effectively_using_the_DelayedArray_framework_for_users.html}{DelayedArray workshop}
#' \href{https://theislab.github.io/zellkonverter/articles/zellkonverter.html}{zellkonverter}
#' @returns Converted single-cell object.
#'
#' @export
#' @examples
#' \dontrun{
#' library(SummarizedExperiment)
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- UpdateSeuratObject(pbmc_small)
#'
#' #### Ingest expression matrix ####
#' sce <- ingest_data(obj=pbmc_small@assays$RNA@counts)
#'
#'#### Ingest EWCElist ####
#' \dontrun{
#' ### Example requires  data from \code{ewceData} package
#' cortex_mrna <- ewceData::cortex_mrna()
#' sce <- ingest_data(obj=cortex_mrna)
#' }
#'
#' #### Ingest Seurat object in memory ####
#' sce <- ingest_data(obj=pbmc_small)
#'
#' #### Ingest HDF5 SingleCellExperiment ####
#' sce <- HDF5Array::saveHDF5SummarizedExperiment(sce, dir = "~/Desktop/pbmc_small_h5", replace=T)
#' ## Read in the sce object directly
#' sce <- ingest_data(obj=sce)
#' ## Read it from disk
#' sce <- ingest_data(obj="~/Desktop/pbmc_small_h5")
#'
#' #### Ingest AnnData ####
#' library(anndata)
#' ## Can point to where anndata is installed (or should be installed)
#' ## Can also just run anndata::install_anndata() and will install via miniconda
#' conda_dir <- dirname(dirname(reticulate::conda_list()[1,]$python))
#' reticulate::use_condaenv(condaenv = conda_dir)
#' reticulate::conda_install(conda = conda_dir, packages = "loompy", pip = T)
#' anndata::install_anndata(method = "conda", conda=conda_dir)
#'
#' # Convert Seurat object to AnnData for example data
#' adata <- anndata::AnnData(X = t(GetAssay(pbmc_small)@counts), obs = pbmc_small@meta.data, var = GetAssay(pbmc_small)@meta.features )
#' ## In memory
#' sce <- ingest_data(obj=adata)
#' ## On disk
#' adata$write_h5ad(filename = "Desktop/pbmc_small.h5ad")
#' sce <- ingest_data(obj = "Desktop/pbmc_small.h5ad")
#'
#'
#' #### Ingest H5Seurat ####
#' library(SeuratDisk)
#' SaveH5Seurat(pbmc_small, filename = "~/Desktop/pbmc_small.h5Seurat", overwrite = T)
#' sce <- ingest_data(obj="~/Desktop/pbmc_small.h5Seurat")
#'
#' #### Ingest loom (from loomR) ####
#' library(loomR)
#' loom <- loomR::create(data=adata, filename = "~/Desktop/pbmc_small.loom", overwrite = T)
#' ## In memory
#' sce <- ingest_data(obj=loom)
#' ## From disk
#' sce <- ingest_data(obj="~/Desktop/pbmc_small.loom")
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


    #### Select output type ####
    output_class <- output_dict(output_class = output_class)
    #### Read ####
    obj <- read_data(path=obj,
                     filetype=input_type,
                     custom_reader=custom_reader,
                     verbose=verbose,
                     ...)
    #### Convert ####
    #### to SummarizedExperiment / SingleCellExperiment ####
    if(is_filetype(output_class,"se")){
        obj_out <- to_se(obj = obj,
                         as_sce = output_class=="singlecellexperiment",
                         verbose = verbose)
    #### to Seurat ####
    } else if(is_filetype(output_class,"seurat") ||
              is_filetype(output_class,"h5seurat")){
        obj_out <- to_seurat(obj = obj,
                             as_h5seurat = is_filetype(output_class,"h5seurat"),
                             verbose = verbose)
    #### to Loom ####
    } else if(is_filetype(output_class,"loom")){
      obj_out <- to_loom(obj = obj,
                         save_path = save_path,
                         verbose = verbose)
    #### to AnnData ####
    }else if(is_filetype(output_class,"anndata")){
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



