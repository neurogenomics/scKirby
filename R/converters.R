#' Converters
#'
#' List all input/output conversions currently supported by \pkg{scKirby}.
#' @export
#' @import data.table
#' @examples
#' map_dt <- converters()
converters <- function(){
  func <- NULL;

  ns <- getNamespaceExports("scKirby")
  to_ns <- grep("_to_",ns,value = TRUE)

  map_dt <- data.table::data.table(func=to_ns)
  map_dt[,c("input","output"):=data.table::tstrsplit(func,"_to_")]
  data.table::setkeyv(map_dt, c("input","output"))
  return(map_dt)
}
