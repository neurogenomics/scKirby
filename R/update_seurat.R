#' Update Seurat object
#'
#' Update a Seurat from V1 --> V3 or V2 --> V3.
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- update_seurat(obj)
update_seurat <- function(obj,
                          verbose=TRUE){

  #### Convert from v1 --> v2 ####
  ## v1 used the class "seurat" all lowercase,
  ## instead of "Seurat" which became the object class name
  ## from Seurat v2 onwards.
  if (methods::is(obj,"seurat")) {
    obj <- seurat1_to_list(obj = obj,
                           verbose = verbose)

  }
  obj <- SeuratObject::UpdateSeuratObject(obj)
  return(obj)
}
