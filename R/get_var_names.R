#' Get var names
#'
#' Get the names of each feature in a single-cell object.
#' @inheritParams converters
#' @returns Vector of observation names.
#'
#' @export
#' @examples
#' obj <- example_obj()
#' vnames <- get_obs_names(obj)
get_var_names <- function(obj) {
  rownames(get_var(obj))
}
