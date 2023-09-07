#' load RData
#'
#' Load processed data using a function that assigns it
#' to a specific variable
#' (so you don't have to guess what the loaded variable name is).
#' @inheritParams base::load
#' @export
load_rdata <- function(file){
  load(file = file)
  get(ls()[ls() != "fileName"])
}
