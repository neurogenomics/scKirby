#' Quantize vector
#'
#' Transform a vector into n quantiles.
#' @param vec A vector.
#' @param n Number of quantiles.
#' @param default Default value to assign to values that fall outside
#' the range of quantiles.
#' @returns A vector of quantiles.
#'
#' @export
#' @examples
#' vec <- runif(100)
#' qvec <- quantize_vector(vec)
quantize_vector <- function(vec,
                            n = 40,
                            default = as.integer(n / 2)
                            ) {

    quantileValues <- rep(0, length(vec))
    breaks <- unique(stats::quantile(vec[vec > 0],
                                     probs = seq(0, 1, by = 1 / n),
                                     na.rm = TRUE
    ))
    if (length(breaks) > 1) {
        quantileValues[vec > 0] <- as.numeric(cut(vec[vec > 0],
                                                  breaks = breaks,
                                                  include.lowest = TRUE
        ))
    } else {
        ## In situations where there's only one non-zero quantile,
        ##  cut() throws an error.
        ## Avoid these situations by using a default quantile.
        messager(
            "+ <2 non-zero quantile bins detected in column.",
            "Assigning these values to default quantile ",
            "(", default, ")"
        )
        quantileValues[vec > 0] <- default
    }
    return(quantileValues)
}
