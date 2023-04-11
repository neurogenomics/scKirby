is_ctd <- function(obj){
  tryCatch({
    utils::getFromNamespace("is_celltypedataset", "EWCE")(ctd = obj)
  }, error=function(e){FALSE})
}
