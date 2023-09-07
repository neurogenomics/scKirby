#' Convert: to \code{DelayedArray}
#'
#' Convert any object to \link[DelayedArray]{DelayedArray} format.
#' @inheritParams converters
#' @inheritParams get_x
#' @inheritDotParams get_x
#' @returns A \link[DelayedArray]{DelayedArray}.
#'
#' @export
#' @importFrom DelayedArray DelayedArray
#' @examples
#' obj <- data.frame(gene=paste0("gene",seq_len(10)),
#'                   matrix(data = 0,nrow=10,ncol = 20))
#' obj2 <- to_delayedarray(obj)
to_delayedarray <- function(obj,
                            as_sparse = TRUE,
                            verbose = TRUE,
                            ...){
  get_x(obj = obj,
        n = 1,
        as_sparse = as_sparse,
        transpose = transpose,
        verbose = verbose,
        ...) |>
    DelayedArray::DelayedArray()
}
