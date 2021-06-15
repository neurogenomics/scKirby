
message_parallel <- function(...){
    system(sprintf('echo "%s"', paste0(..., collapse="")))
}

messager <- function(..., v=T){if(v){message(paste(...))}}

printer <- function(..., v=T){if(v){print(paste(...))}}


#' loadRData
#'
#' Load processed data using a function that assigns it
#' to a specific variable (so you don't have to guess what the loaded variable name is).
#'
#' @export
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}




assign_cores <- function(worker_cores=.90,
                         verbose=T){
    # Enable parallelization of HDF5 functions
    ## Allocate ~10% of your available cores to non-parallelized processes
    total_cores <- parallel::detectCores()
    if(worker_cores<1){
        reserved_cores <-  ceiling(total_cores*(1-worker_cores))
        workers <- total_cores - reserved_cores
    } else {
        workers <- worker_cores
        reserved_cores <-  total_cores - workers
    }
    message("+ ",workers," core(s) assigned as workers (",reserved_cores," reserved).")
    DelayedArray::setAutoBPPARAM(BiocParallel::MulticoreParam(workers))
    DelayedArray:::set_verbose_block_processing(verbose)
    return(list(worker_cores=workers,
                reserved_cores=reserved_cores,
                total_cores=total_cores))
}


