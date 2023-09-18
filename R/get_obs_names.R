#' Get obs names
#'
#' Get the names of each observation in a single-cell object.
#' @inheritParams converters
#' @returns Vector of observation names.
#'
#' @export
#' @examples
#' obj <- example_obj()
#' obs_names <- get_obs_names(obj)
get_obs_names <- function(obj) {
  rownames(get_obs(obj))
}
