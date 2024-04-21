#' Converters
#'
#' Functions to convert one single-cell format to another.
#' @param obj A single-cell object supported by \pkg{scKirby}.
#'  See \link{converters} for a table of all supported conversions.
#' @param verbose Print messages.
#' @param ... Additional arguments passed to the respective converter function.
#' \link{converters} function.
#'
#' @export
#' @examples
#' map_dt <- converters()
converters <- function(obj, verbose, ...){
  func <- NULL;

  ns <- getNamespaceExports("scKirby")
  to_ns <- grep("_to_",ns,value = TRUE)
  map_dt <- data.table::data.table(func=to_ns)
  map_dt[,c("input","output"):=data.table::tstrsplit(func,"_to_")]
  data.table::setkeyv(map_dt, c("input","output"))
  return(map_dt)
}
