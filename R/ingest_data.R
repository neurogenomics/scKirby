#' Import and standardize scRNAseq data across different formats
#'
#' Automatically infers data format of scRNAseq object, or a path to that object.
#' It then uses the appropriate functions to import that data and convert it to a
#' \pkg{SingleCellExperiment}, which is recognized by other \pkg{EWCE} functions.
#'
#' @param obj Single-cell data object or path to saved single-cell data.
#' @param input_type Format of \code{obj}. By default, the type will be inferred.
#' @param output_type Format to convert \code{obj} to.
#' @param custom_reader Custom function to read \code{obj} into R.
#' @param save_dir Directory to save the converted \code{obj}.
#' @param filename Name to save the converted \code{obj}.
#' @param save_output Whether or not to save the converted \code{obj}.
#' @param overwrite If a file of the same name exists, overwrite it.
#' @param return_filepath If \code{TRUE}, a list with both the converted object
#'  and the saved file path will be returned (instead of just the converted object).
#' @param ... Additional arguments to be passed to \code{scKirby::read_data()}.
#'
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
#' @import dplyr
#' @source
#' \href{https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html}{SeuratDisk}
#' \href{https://github.com/rcannood/anndata}{anndata (R)}
#' \href{https://anndata.readthedocs.io/en/latest/}{anndata (python)}
#' \href{https://satijalab.org/loomR/loomR_tutorial.html}{loomR}
#' \href{https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html}{SingleCellExperiment}
#' \href{https://petehaitch.github.io/BioC2020_DelayedArray_workshop/articles/Effectively_using_the_DelayedArray_framework_for_users.html}{DelayedArray workshop}
#' \href{https://theislab.github.io/zellkonverter/articles/zellkonverter.html}{zellkonverter}
#' @export
ingest_data <- function(obj,
                        input_type="guess",
                        output_type="SingleCellExperiment",
                        custom_reader=NULL,
                        save_dir=tempdir(),
                        filename="scKirby_output",
                        save_output=T,
                        overwrite=F,
                        return_filepath=F,
                        verbose=T,
                        ...){
    # output_type="SingleCellExperiment";custom_reader=NULL;save_dir=tempdir();filename="scKirby_output"; overwrite=F; return_filepath=F;verbose=T;
    cdict <- class_dictionary()
    output_types <- list(sce=tolower(cdict$sce),
                         hdf5se=tolower(cdict$hdf5se),
                         seurat=tolower(cdict$seurat),
                         h5seurat=tolower(cdict$h5seurat)
                         # loom=c("loom"),
                         # anndata=c("anndata","h5ad")
                         )
    if(!tolower(output_type) %in% unname(unlist(output_types))){
        stop("output_type must be one of the following: ", paste(unname(unlist(output_types)), collapse = ", "))
    }

    #### Read ####
    # Separate the reading/conversion process
    ## bc you don't always know what kind of data you're reading in (esp .rds/.rda files).
    object <- read_data(obj=obj,
                         filetype=input_type,
                         custom_reader=custom_reader,
                         save_dir=save_dir,
                         overwrite=overwrite,
                         verbose=verbose,
                         ...)

    #### Convert ####
    ##  to SingleCellExperiment
    if(tolower(output_type) %in% output_types$sce){
        object_out <-  to_sce(object = object,
                              verbose = verbose)
        if(tolower(output_type)==tolower("SummarizedExperiment")){
            object_out <- sce_to_se(object=object_out,
                                    verbose=verbose)
        }
    }
    ## to Seurat
    if(tolower(output_type) %in% c(output_types$seurat, output_types$h5seurat)){
        object_out <- to_seurat(object=object,
                                save_dir=save_dir,
                                verbose=verbose)
    }

    #### Save ####
    if(save_output){
        filepath <- save_data(object=object_out,
                              output_type=output_type,
                              save_dir=save_dir,
                              filename=filename,
                              overwrite=overwrite,
                              verbose=verbose)
    } else {filepath <- NULL }

    #### Return ####
    if(return_filepath) return(filepath=filepath, object=object_out) else return(object_out)
}



