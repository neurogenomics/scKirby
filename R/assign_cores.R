assign_cores <- function(workers=.90,
                         return_list=FALSE,
                         verbose=TRUE){

  # Enable parallelization of HDF5 functions
  ## Allocate ~10% of your available cores to non-parallelized processes
  total_cores <- parallel::detectCores()
  if(workers<1){
    reserved_cores <-  ceiling(total_cores*(1-workers))
    workers <- total_cores - reserved_cores
  } else {
    workers <- workers
    reserved_cores <-  total_cores - workers
  }
  message("+ ",workers," core(s) assigned as workers (",reserved_cores," reserved).")
  DelayedArray::setAutoBPPARAM(BiocParallel::MulticoreParam(workers))
  DelayedArray:::set_verbose_block_processing(verbose)
  if(isTRUE(return_list)){
    return(list(workers=workers,
                reserved_cores=reserved_cores,
                total_cores=total_cores))
  }
}
