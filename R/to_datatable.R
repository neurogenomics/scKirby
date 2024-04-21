#' To \code{data.table}
#'
#' @describeIn converters
#' Convert any single-cell object to
#' \link[data.table]{data.table} format.
#' @inheritParams converters
#' @inheritParams get_x
#' @returns A \code{data.table} object.
#'
#' @export
#' @examples
#' obj <- example_obj("cds")
#' obj2 <- to_dataframe(obj)
to_datatable <- function(obj,
                         n = 1,
                         verbose = TRUE,
                         ...){

  #### Check if class is supported ####
  check_supported(obj)
  ata.table::as.data.table(
    x = get_x(obj = obj,
              n = n,
              verbose = verbose,
              ...),
    keep.rownames = TRUE)
}
