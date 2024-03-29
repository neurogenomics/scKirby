dict_suffix <- function(){
  list(rds=".rds",
       rdata=c(".rda",".rdata"),
       matrix=".mtx",
       data.table=c(".csv",".tsv",".csv.gz",".tsv.gz",".txt",".txt.gz"),
       anndata=".h5ad",
       seurat=".rds",
       h5seurat=".h5seurat",
       loom=".loom",
       h5=c(".h5$",".sce"))
}
