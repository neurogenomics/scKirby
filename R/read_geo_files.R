
read_geo_files <- function(file_dir,
                           verbose=TRUE){
  # file_dir <- "~/projects/model_celltype_conservation/raw_data/scRNAseq/Raj2020/"
  file_prefixes <- list.files(file_dir,
                              pattern = "-barcodes.tsv.gz") |>
    gsub(pattern = "-barcodes.tsv.gz",replacement = "")

  core_allocation <- assign_cores(workers = .9, return_list = TRUE)

  sce_list <- parallel::mclapply(file_prefixes, function(x){
    messager(x, parallel = TRUE, v = verbose)
    #### Read in components ####
    colDat <- read.delim2(file.path(file_dir,paste0(x,"-barcodes.tsv.gz")),
                          header = FALSE,
                          col.names = "barcodes")
    colDat$sample <- x
    rowDat <- read.delim2(file.path(file_dir,paste0(x,"-genes.tsv.gz")),
                          header = FALSE,
                          col.names = c("ensembl_id","gene_symbol"))
    mat <- DelayedArray::DelayedArray(
      Matrix::readMM(file.path(file_dir,paste0(x,"-matrix.mtx.gz")))
    )
    #### Construct SCE object ####
    sce_sub <- SingleCellExperiment::SingleCellExperiment(
      assays = list(raw = mat),
      colData = colDat,
      rowData = rowDat
    )
    sce_sub <- check_se_rownames(sce_sub,
                                  rownames_var = "gene_symbol",
                                  verbose = FALSE)
    return(sce_sub)
  }, mc.cores = core_allocation$workers)
  messager("Merging all SCE objects into one.",v=verbose)
  sce <- do.call("cbind",sce_list)
  return(sce)
}



