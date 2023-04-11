#' load RData
#'
#' Load processed data using a function that assigns it
#' to a specific variable (so you don't have to guess what the loaded variable name is).
#'
#' @export
load_rdata <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
