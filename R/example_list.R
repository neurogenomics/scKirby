example_list <- function(obj=SeuratObject::pbmc_small,
                         as_paths=FALSE,
                         verbose=TRUE){

  obj <- to_list(obj = obj,
                 verbose = verbose)
  if(isTRUE(as_paths)){
    #### Save data ####
    obj <- lapply(stats::setNames(names(obj),
                                   names(obj)),
                   function(nm){
      tmp <- tempfile(fileext = paste0("_example_",nm,".rds"))
      saveRDS(obj[[nm]],tmp)
      tmp
    })
  }
  return(obj)
}
