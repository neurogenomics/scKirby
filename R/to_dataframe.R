#' To \code{data.frame}
#'
#' @describeIn converters
#' Convert any single-cell object to
#' \link[base]{data.frame} format.
#' @inheritParams converters
#' @inheritParams get_x
#' @returns A \code{data.frame} object.
#'
#' @export
#' @examples
#' obj <- example_obj("cds")
#' obj2 <- to_dataframe(obj)
to_dataframe <- function(obj,
                         n = 1,
                         verbose = TRUE,
                         ...){

  #### Check if class is supported ####
  check_supported(obj)
  as.data.frame(
    x = get_x(obj = obj,
              n = n,
              verbose = verbose,
              ...),
    keep.rownames = TRUE)
}
