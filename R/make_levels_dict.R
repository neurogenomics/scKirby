make_levels_dict <- function(obs,
                             level_cols){
  #### Uns dict ####
  uns = data.frame(level = names(level_cols),
                   label = sapply(level_cols,paste,collapse="."))
  #### Get number of unique values per level ####
  uns["n"] = sapply(uns$level, function(x){length(unique(obs[[x]]))})
  return(uns)
}
