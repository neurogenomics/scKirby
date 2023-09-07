get_obsm_mofa <- function(obj,
                          verbose=TRUE){
  #### Add MOFA factors ####
  obsm <- obj@expectations$Z
  if(is.list(obsm) && length(obsm)>1){
    messager(">1 group found in MOFA embeddings.",
             "Using only the first:",shQuote(names(obsm)[1]),v=verbose)
    obsm <- obsm[[1]]
    if(!is.null(obj@samples_metadata$sample)){
      rownames(obsm) <- obj@samples_metadata$sample
    }
  }
  obsm <- list(mofa=as.matrix(as.data.frame(obsm)))
  #### Add further reduced factors ####
  if(!is.null(obj@dim_red)){
    for(nm in names(obj@dim_red)){
      nm2 <- paste("mofa",tolower(nm),sep="_")
      obsm[[nm2]] <- as.matrix(obj@dim_red[[nm]][-1])
      if(!is.null(obj@samples_metadata$sample)){
        rownames(obsm[[nm2]]) <- obj@samples_metadata$sample
      }
    }
  }
  return(obsm)
}
