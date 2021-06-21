save_sce <- function(object,
                     save_dir=tempdir(),
                     filename="scKirby_output",
                     verbose=T){
  filename <- file.path(save_dir, paste0(filename,".rds"))
  messager("+ Saving SingleCellExperiment:", filepath,v=verbose)
  dir.create(path = dirname(filename), showWarnings = F, recursive = T)
  saveRDS(object, file = filename)
  return(filename)
}
