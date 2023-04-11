is_filetype <- function(x,
                        group=NULL){
  if(is.null(x)) return(FALSE)
  dict <- filetype_dict()
  matches <- lapply(dict, function(y){
    grepl(paste(paste0(gsub("\\.","\\.",y),"$"),collapse = "|"),x)
  })
  groups <- names(matches[unlist(matches)])
  if(is.null(group)){
    return(groups)
  } else {
    group <- tolower(group)
    if(!group %in% names(dict)){
      stopper(group,"is not a name in filetype_dict.")
    }
    return(group %in% groups)
  }
}
