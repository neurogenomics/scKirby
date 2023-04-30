check_keys <- function(all_keys,
                            keys){
  if(length(keys)>0){
    all_keys <- all_keys[tolower(all_keys) %in% tolower(keys)]
  }
  return(all_keys)
}
