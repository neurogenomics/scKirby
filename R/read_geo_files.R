
read_geo_files <- function(file_dir){
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



