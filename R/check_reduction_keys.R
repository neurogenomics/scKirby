check_reduction_keys <- function(all_keys,
                                 keys){
  if(!is.null(keys)){
    all_keys <- all_keys[tolower(all_keys) %in% tolower(keys)]
  }
  return(all_keys)
}
