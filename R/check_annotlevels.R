check_annotlevels <- function(obj,
                              annotLevels){

  if(!is.list(annotLevels)){
    stopper("annotLevels must be a list.")
  }
  obs <- get_obs(obj)
  dims <- dim(obj)
  if(is_class(obj,"anndata")) dims <- rev(dims)
  #### Check annotLevels ####
  if(all(mapply(is.character,annotLevels))){
    # Check
    check <- lapply(seq_len(length(annotLevels)), function(i){
      if(length(annotLevels[[i]])!=dims[2]){
        stopper("Element",i,"in annotLevels must have the same length",
                "as the number of cells in obj.")
      }
    })
  } else if(all(names(annotLevels) %in% names(obs))){
    messager("Extracting annotLevels from obs")
    annotLevels <- lapply(stats::setNames(names(annotLevels),
                                          names(annotLevels)),
           function(x){
           as.character(obs[[x]])
           })
  }
  #### Check for blanks ####
  annotLevels <- lapply(annotLevels, function(x){
    if("" %in% unique(x)){
      x[x==""] <- "NA"
    }
    return(x)
  })
  return(annotLevels)
}
