#' Quantize a data.table
#'
#' @param dat A data.table
#' @param vars A character vector of column names
#' @param n The number of quantiles to use
#' @returns A data.table with the same columns as \code{dat} but with
#'  the values quantized into \code{n} quantiles.
#'
#' @export
#' @examples
#' dat <- data.table::data.table(x = rnorm(100))
#' quantize_dt(dat, vars = "x")
#' print(dat)
quantize_dt <- function(dat,
                        vars,
                        n=40){

    if(methods::is(dat,"data.frame")&&
       !methods::is(dat,"data.table")){
        dat <- data.table::as.data.table(dat)
    }
    for(v in vars){
        dat[,(v):=quantize_vector(get(v),n=n)]
    }
}
