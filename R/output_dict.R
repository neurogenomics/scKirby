#' Output dictionary
#'
#' A list of all object classes that the single-cell dataset can be returned as.
#' Classes are grouped by families of related classes.
#' @param output_class Name of the output class you wish to produce.
#' Alternatively, set to \code{NULL} to return all possible output classes.
#' @returns A valid \code{output_class}.
#'
#' @export
#' @examples
#' output_class <- output_dict(output_class="seurat")
output_dict <- function(output_class=NULL){

  cdict <- class_dict()
  output_classes <- list(se=tolower(cdict$se),
                       hdf5se=tolower(cdict$hdf5se),
                       seurat=tolower(cdict$seurat),
                       h5seurat=tolower(cdict$h5seurat),
                       loom=c("loom"),
                       list="list"
                       # anndata=c("anndata","h5ad")
                       )
  if(is.null(output_class)){
   return(output_classes)
  } else {
    output_class <- tolower(output_class)[1]
    if(!output_class %in% unname(unlist(output_classes))){
      stopper("output_type must be one of the following:",
           paste(unname(unlist(output_classes)), collapse = ", "))
    } else {
      return(output_class)
    }
  }
}
