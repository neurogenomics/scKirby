


#' Import and standardize scRNAseq data across different formats
#'
#' Automatically infers data format of scRNAseq object, or a path to that object.
#' It then uses the appropriate functions to import that data and convert it to a
#' \pkg{SingleCellExperiment}, which is recognized by other \pkg{EWCE} functions.
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
#' \href{https://github.com/cellgeni/sceasy}{sceasy}
#' \href{https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html}{SeuratDisk}
#' \href{https://github.com/rcannood/anndata}{anndata (R)}
#' \href{https://anndata.readthedocs.io/en/latest/}{anndata (python)}
#' \href{https://satijalab.org/loomR/loomR_tutorial.html}{loomR}
#' \href{https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html}{SingleCellExperiment}
#' \href{https://petehaitch.github.io/BioC2020_DelayedArray_workshop/articles/Effectively_using_the_DelayedArray_framework_for_users.html}{DelayedArray workshop}
#' \href{https://theislab.github.io/zellkonverter/articles/zellkonverter.html}{zellkonverter}
#' @export
ingest_data <- function(obj,
                        filetype="guess",
                        custom_reader=NULL,
                        sce_save_dir=NULL,
                        quicksave_HDF5=T,
                        overwrite=F,
                        verbose=T,
                        ...){
    # Separate the reading/SCE conversion process
    ## bc you don't always know what kind of data you're reading in
    ### (esp .rds/.rda files).
    object <- read_scRNAseq_data(obj=obj,
                                 filetype=filetype,
                                 custom_reader=custom_reader,
                                 sce_save_dir=sce_save_dir,
                                 overwrite=overwrite,
                                 verbose=verbose,
                                 ...)
    sce <- convert_to_SCE(object = object,
                          verbose = verbose)
    sce <- save_SCE(sce=sce,
                    save_dir=sce_save_dir,
                    quicksave_HDF5=quicksave_HDF5,
                    overwrite=overwrite,
                    verbose=verbose)
    return(sce)
}




read_scRNAseq_data <- function(obj,
                               filetype="guess",
                               custom_reader=NULL,
                               sce_save_dir=NULL,
                               overwrite=F,
                               verbose=T,
                               ...){
    if(!is.null(custom_reader)){
        messager("+ Reading in with custom_reader function...")
        object <- custom_reader(obj, ...)
        return(object)
    }

    if(class(obj)[1]=="character"){
        messager("+ Reading from disk...")
        #### Generic RDS ####
        if(endsWith(tolower(obj), suffix=".rds") | tolower(filetype)=="rds"){
            messager("+ Reading in .rds file of unknown type...")
            object <- readRDS(obj, ...)
            return(object)
        }
        #### Generic RDS ####
        if(any(endsWith(tolower(obj), suffix=c(".rda",".rdata"))) | tolower(filetype)=="rda"){
            messager("+ Reading in .rda file of unknown type...")
            object <- loadRData(obj)
            return(object)
        }
        #### .mtx folder ####
        if((endsWith(tolower(obj), suffix=".mtx") & dir.exists(obj)) | tolower(filetype)=="mtx_dir"){
            messager("+ Matrix folder format (.mtx) detected. Importing as Seurat object...")
            object <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(data.dir=obj, ...))
            return(object)
        }
        #### .mtx matrix ####
        if((endsWith(tolower(obj), suffix=".mtx") & (!dir.exists(obj)) ) | tolower(filetype)=="mtx"){
            messager("+ Expression matrix (.mtx) detected. Importing as sparse dgCMatrix...",v=verbose)
            object <- data.table::fread(obj, stringsAsFactors = F,
                                        data.table = F)
            object <- object %>% tibble::column_to_rownames(colnames(object)[1]) %>%
                as.matrix() %>% Matrix::Matrix(sparse=T)
            return(object)
        }
        #### .csv/.tsv matrix ####
        if((any(endsWith(tolower(obj), suffix=c(".csv",".tsv",".csv.gz",".tsv.gz")))) | tolower(filetype) %in% c("csv","tsv")){
            messager("+ Expression matrix (.csv|.tsv) detected. Importing as sparse dgCMatrix...",v=verbose)
            object <- data.table::fread(obj, stringsAsFactors = F,
                                        data.table = F)
            object <- object %>% tibble::column_to_rownames(colnames(object)[1]) %>%
                as.matrix() %>% Matrix::Matrix(sparse=T)
            return(object)
        }
        #### AnnData ####
        if((endsWith(tolower(obj), suffix=".h5ad")) | tolower(filetype)=="h5ad"){
            messager("+ AnnData format (.h5ad) detected. Importing as AnnData object...",v=verbose)
            #### anndata method
            # Anndata adds another dependency, but at least it works unlike
            object <- anndata::read_h5ad(filename = obj)

            #### sceasy method
            # object <- sceasy::convertFormat(obj, from="anndata", to="sce")

            #### Seurat method
            ## This is now deprecated, and for some reason doesn't offer back compatibility by calling to SeuratDisk...
            # object <- Seurat::ReadH5AD(obj, ...)

            ## SeuratDisk is currently broken, with no word from the developer...
            ## https://github.com/mojaveazure/seurat-disk/issues/41
            # object <- SeuratDisk::Convert(source = obj,
            #                               dest = sce_save_dir,
            #                               overwrite = overwrite,
            #                               verbose = verbose,
            #                               ...)
            return(object)
        }
        #### H5Seurat ####
        if((endsWith(tolower(obj), suffix=".h5seurat")) | tolower(filetype)=="h5seurat"){
            messager("+ h5Seurat format (.h5Seurat) detected. Importing as Seurat object...",v=verbose)
            object <- SeuratDisk::LoadH5Seurat(file = obj, ...)
            return(object)
        }
        #### Loom ####
        if((endsWith(tolower(obj), suffix=".loom")) | tolower(filetype)=="loom"){
            messager("+ Loom format (.loom) detected. Importing as SingleCellLoomExperiment object...",v=verbose)
            #### anndata method
            ## anndata::read_loom has difficulties identifying right loompy location.
            # anndata::read_loom(filename=obj, validate=F, ...)

            #### loomR method
            ## skip.validate must =F, or else you won't be able to extract the matrix
            # object <- loomR::connect(filename=obj, skip.validate = F)

            ### SeuratDisk method
            object <- SeuratDisk::LoadLoom(obj)

            #### sceasy method
            # rhdf5::h5disableFileLocking() ## Causes error otherwise
            # object <- sceasy::convertFormat(obj, from="loom", to="sce", ...)
                                            # outFile = gsub(".loom",".sce.rds",obj))
            return(object)
        }
        #### HDF5Array SummarizedExperiment/SingleCellExperiment ####
        if((any(endsWith(tolower(obj), suffix=c("h5","sce"))) |
            tolower(filetype) %in% c("HDF5Array","SummarizedExperiment","SingleCellExperiment","SingleCellLoomExperiment") ) &
            dir.exists(obj) ){
                if(file.exists(file.path(obj,"assays.h5")) & file.exists(file.path(obj,"se.rds")) ){
                    messager("+ HDF5Array format (.h5) detected. Importing as SingleCellExperiment object...",v=verbose)
                    object <- HDF5Array::loadHDF5SummarizedExperiment(obj, ...)
                    return(object)
                }
        }

        # if( endsWith(tolower(obj), suffix=c(".h5")) & (!dir.exists(obj))){
        #     # When .h5 is a file, not a folder
        #     # obj <- "/Volumes/bms20/projects/neurogenomics-lab/live/GitRepos/model_celltype_conservation/raw_data/scRNAseq/LaManno2020/LaManno2020_sparse.h5"
        #     # rhdf5::h5closeAll()
        #     # h5_contents <- rhdf5::h5ls(obj, datasetinfo=F)
        #     # object <- rhdf5::h5read(file =obj, name = "/")
        #     object <- anndata::read_hdf(obj, key = "/")
        # }

    } else {
        messager("+ Returning object directly...",v=verbose)
        return(obj)
    }
}





read_GEO_files <- function(file_dir){
    # file_dir <- "~/projects/model_celltype_conservation/raw_data/scRNAseq/Raj2020/"
    file_prefixes <- list.files(file_dir,
                                pattern = "-barcodes.tsv.gz") %>%
        gsub(pattern = "-barcodes.tsv.gz",replacement = "")

    core_allocation <- assign_cores(worker_cores = .9)

    sce_list <- parallel::mclapply(file_prefixes, function(x){
        message_parallel(x)
        #### Read in components ####
        colDat <- read.delim2(file.path(file_dir,paste0(x,"-barcodes.tsv.gz")),
                              header = F, col.names = "barcodes")
        colDat$sample <- x
        rowDat <- read.delim2(file.path(file_dir,paste0(x,"-genes.tsv.gz")),
                              header = F, col.names = c("ensembl_id","gene_symbol"))
        mat <- DelayedArray::DelayedArray(Matrix::readMM(file.path(file_dir,paste0(x,"-matrix.mtx.gz"))))
        #### Construct SCE object ####
        sce_sub <- SingleCellExperiment::SingleCellExperiment(
            assays      = list(raw = mat),
            colData     = colDat,
            rowData     = rowDat
        )
        sce_sub <- check_sce_rownames(sce_sub,
                                      rownames_var = "gene_symbol",
                                      verbose = F)
        return(sce_sub)
    }, mc.cores = core_allocation$worker_cores)
    printer("Merging all SCE objects into one...")
    sce <- do.call("cbind",sce_list)
    return(sce)
}





convert_to_SCE <- function(object,
                           verbose=T){
    messager("Converting formats:",v=verbose)
    core_allocation <- assign_cores(worker_cores = .90)

    ewce_classes <- c("EWCE_list","SCE_list")
    matrix_classes <- c("data.table","data.frame","tbl_df","tbl","matrix","Matrix","array","DelayedArray","DelayedMatrix",names(getClass("Matrix")@subclasses)) # more than 40 ..
    loom_classes <-  c("loom","H5File","H5RefClass")
    sce_classes <- c("SingleCellLoomExperiment","SingleCellExperiment","SummarizedExperiment")
    anndata_classes <- c("AnnData","AnnDataR6","AnnDataR6R6")
    seurat_classes <- c("Seurat")
    supported_classes <- c(matrix_classes, loom_classes, sce_classes, anndata_classes, seurat_classes, ewce_classes)
    supported_classes_print <- c("matrix (any subclass)",loom_classes, sce_classes, anndata_classes, seurat_classes, ewce_classes)

    #### EWCE_list class ####
    if(class(object)[1]=="list" & all(c("exp","annot") %in% names(object)) ){
        class(object) <- "EWCE_list"
    }
    #### CTD class ####
    if(class(object)[1]=="list" & length(object)>1){
        if(all(c("mean_exp","specificity") %in% names(object[[1]]) )){
            messager("+ CTD ==> SingleCellExperiment",v=verbose)
            # object <- readRDS("~/Desktop/model_celltype_conservation/processed_data/EWCE/CTD_list.rds")[[1]]
            ctd <- object
            #### Name CTD levels ####
            if(is.null(names(ctd))){
                names(ctd) <- paste0("level_",1:length(ctd))
            } else {
                names(ctd) <- names(ctd)
            }
            sce_list <- lapply(names(ctd), function(lvl){
                message("Converting level: ",lvl)
                ctd_lvl <- ctd[[lvl]]
                #### Use matrices that are present ###
                matrix_list <- list()
                for(mtx_name in c("mean_exp","median_exp",
                                  "specificity","median_specificity","specificity_quantiles")){
                    if(mtx_name %in% names(ctd_lvl)){ matrix_list[[mtx_name]] <- DelayedArray::DelayedArray(as(as(ctd_lvl[[mtx_name]], "matrix"), "sparseMatrix")) }
                }
                sce <- SingleCellExperiment::SingleCellExperiment(
                    assays      =  matrix_list,
                    colData     =  data.frame(colnames(matrix_list[[1]])) %>% `colnames<-`(lvl),
                    rowData     =  data.frame(gene=row.names(matrix_list[[1]]), row.names = row.names(matrix_list[[1]]))
                )
                sce <- check_sce_rownames(sce, verbose = verbose)
            }) %>% `names<-`(names(ctd))
            ## "SCE_list" class messes up other functions that expect class "list"
            # class(sce_list) <- "SCE_list"
            return(sce_list)
        }
    }




    #### Check if class is supported ####
    if(!class(object)[1] %in% supported_classes){
        stop("Unsupported class detected: ",class(object),"\n\n",
             "Data object must be at least one of the following classes:\n\n",
             paste(supported_classes_print, collapse = ", "))
    }
    ### EWCE_list ####
    if(class(object)[1] == "EWCE_list"){
        messager("+ EWCE_list ==> SingleCellExperiment",v=verbose)
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays      = list(raw = DelayedArray::DelayedArray(as(as.matrix(object$exp), "sparseMatrix"))),
            colData     =  object$annot,
            rowData     =  row.names(object$exp)
        )
        sce <- check_sce_rownames(sce, verbose = verbose)
        return(sce)
    }

    #### Matrices ####
    if(class(object)[1] %in% matrix_classes){
        messager("+ Matrix ==> SingleCellExperiment",v=verbose)
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays      = list(raw = DelayedArray::DelayedArray(as(as.matrix(object), "sparseMatrix"))),
        )
        sce <- check_sce_rownames(sce, verbose = verbose)
        return(sce)
    }
    #### Seurat ####
    if(class(object)[1] %in% seurat_classes){
        messager("+ Seurat ==> SingleCellExperiment",v=verbose)
        # object <- Seurat::pbmc_small ## example
        ### Seurat::as.SingleCellExperiment() is currently broken
        # sce <- Seurat::as.SingleCellExperiment(object)

        ### Must convert manually instead
        assay1 <- names(object@assays)
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays      = list(raw = DelayedArray::DelayedArray(as(object@assays[[assay1]]@counts, "sparseMatrix"))),
            colData     = object@meta.data,
            rowData     = object@assays[[assay1]]@meta.features
        )
        sce <- check_sce_rownames(sce, verbose = verbose)
        return(sce)
    }
    #### AnnData ####
    if(class(object)[1] %in% anndata_classes){
        messager("+ AnnData ==> SingleCellExperiment",v=verbose)
        # sceasy::convertFormat(obj = object, from = "anndata", to="seurat",
        #                       outFile = "~/Desktop/tmp.rds")
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays      = list(raw = DelayedArray::DelayedArray(as( Matrix::t(object$X), "sparseMatrix"))),
            colData     = object$obs,
            rowData     = object$var
        )
        sce <- check_sce_rownames(sce, verbose = verbose)
        return(sce)
    }
    if(class(object)[1] %in% loom_classes){
        messager("+ loom ==> SingleCellExperiment",v=verbose)
        # Import as a Seurat object first for convenience
        object <- SeuratDisk::LoadLoom(file = object$filename)
        assay1 <- names(object@assays)
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays      = list(raw = DelayedArray::DelayedArray(as(object@assays[[assay1]]@counts, "sparseMatrix"))),
            colData     = object@meta.data,
            rowData     = object@assays[[assay1]]@meta.features
        )
        sce <- check_sce_rownames(sce, verbose = verbose)
        return(sce)
    }
    #### SingleCellExperiment/SummarizedExperiment ####
    if(class(object)[1] %in% sce_classes){
        messager("+ == SummarizedExperiment",v=verbose)
        sce <- object
        sce <- check_sce_rownames(sce, verbose = verbose)
        return(sce)
    }
}


check_sce_rownames <- function(sce,
                               rownames_var=NULL,
                               remove_duplicates=T,
                               verbose=T){
    printer("Checking SCE rownames...",v=verbose)
    rowDat <- SummarizedExperiment::rowData(sce)
    rownames_var <- if(is.null(rownames_var)) "Gene" else rownames_var
    if( is.null( S4Vectors::rownames(rowDat)) ){
        if(rownames_var %in% colnames(rowDat)){
            printer("+ Assigning rownames to:",rownames_var,v=verbose)
            # IMPORTANT! only S4Vectors::rownames can assign row names.
            S4Vectors::rownames(sce) <- rowDat[[rownames_var]] # Automatically assigns same rownames to assay and rowData
        }else {
            warning("Cannot identify rownames. Please set rownames (rowDat) first with: S4Vectors::rownames(sce) <- gene_names")
        }
    }
    if(remove_duplicates){
        printer("+ Removing duplicate gene rows.",v=verbose)
        sce <- sce[!base::duplicated(S4Vectors::rownames(sce)),]
    }
    return(sce)
}


check_sce_colnames <- function(sce,
                               colnames_var=NULL,
                               verbose=T){
    printer("Checking SCE colnames.",v=verbose)
    colDat <- SummarizedExperiment::colData(sce)
    if(!is.null(colnames_var)){
        if(colnames_var %in% colnames(colDat)){
            printer("Assigning colnames to:", colnames_var,v=verbose)
            S4Vectors::colnames(sce) <- colDat[[colnames_var]]
        } else {
            printer(colnames_var,"not found in colData.",v=verbose)
        }
    }
    return(sce)
}


save_SCE <- function(sce,
                     save_dir,
                     quicksave_HDF5=T,
                     overwrite=F,
                     verbose=T){
    if(!is.null(save_dir)){
        if(quicksave_HDF5 & file.exists(file.path(save_dir,"assays.h5")) & (overwrite==F)){
            messager("+ Updating existing HDF5...",v=verbose)
            sce <- tryCatch(expr = {
                        HDF5Array::quickResaveHDF5SummarizedExperiment(x=sce, verbose=verbose)
                    },
                    error=function(e){
                        HDF5_hard_save(sce=sce,
                                       save_dir=save_dir,
                                       overwrite=overwrite,
                                       verbose=verbose)
            })
        } else {
            if( (!dir.exists(save_dir)) | (overwrite) ){
                sce <- HDF5_hard_save(sce=sce,
                                      save_dir=save_dir,
                                      overwrite=overwrite,
                                      verbose=verbose)

            } else {
                messager("+ Returning existing SCE object.",v=verbose)
            }
        }
    }
    return(sce)
}


HDF5_hard_save <- function(sce,
                           save_dir,
                           overwrite=T,
                           verbose=T){
    messager("+ Writing new HDF5...",v=verbose)
    # DON'T create the HDF5 dir itself (will return an error about overwriting)
    dir.create(dirname(save_dir), showWarnings = F, recursive = T)
    # IMPORTANT!: set as.sparse=T if you have the latest version of HDF5Array (1.8.11)
    pkg_ver <- packageVersion("HDF5Array")
    pkg_ver_split <- strsplit(as.character(pkg_ver),".", fixed = T)[[1]]
    pkg_V <- as.numeric(paste(pkg_ver_split[1], pkg_ver_split[2], sep="."))
    if(pkg_V>=1.8){
        sce <- HDF5Array::saveHDF5SummarizedExperiment(x=sce,
                                                       dir=save_dir,
                                                       verbose=verbose,
                                                       as.sparse=T,
                                                       replace=overwrite)
    }else{
        sce <- HDF5Array::saveHDF5SummarizedExperiment(x=sce,
                                                       dir=save_dir,
                                                       verbose=verbose,
                                                       replace=overwrite)
    }
    return(sce)
}



