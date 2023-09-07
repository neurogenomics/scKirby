#' Convert: \link[base]{list} ==> \code{DimReducObject}
#'
#' Create a \pkg{Seurat} \code{DimReducObject}
#' from embeddings and loadings.
#' @param obj A named list containing:
#' \itemize{
#' \item{obsm : }{Sample embeddings}
#' \item{varm : }{Feature loadings}
#' \item{stdev : Standard deviation (if applicable) for the
#'  dimensional reduction}
#' }
#' @param assay Assay to use.
#' @param key Key to use (name of embedding).
#' @inheritDotParams SeuratObject::CreateDimReducObject
#'
#' @export
#' @importFrom SeuratObject CreateDimReducObject
list_to_dimreduc <- function(obj,
                             assay,
                             key,
                             ...) {

    Seurat::CreateDimReducObject(
        embeddings = obj$obsm,
        loadings = obj$varm,
        stdev = if (is.null(obj$stdev)) {
            numeric()
        } else {
            as.numeric(obj$stdev)
        },
        assay = assay,
        key = key,
        ...)
}
