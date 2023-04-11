is_suffix <- function(x,
                      group=NULL){
  dict <- suffix_dict()
  matches <- lapply(dict, function(y){
    grepl(paste(paste0(y,"$"),collapse = "|"),x)
  })
  groups <- names(matches[unlist(matches)])
  if(is.null(group)){
    return(groups)
  } else {
    group <- tolower(group)
    if(!group %in% names(dict)){
      stopper(group,"is not a name in suffix_dict.")
    }
    return(group %in% groups)
  }
}
