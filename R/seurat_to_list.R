#' Convert: \code{Seurat} ==> \code{list}
#'
#' @export
#' @importFrom methods .hasSlot slot
#' @importFrom stats setNames
#' @examples
#' obj <- example_obj("seurat")
#' obj2 <- seurat_to_list(obj)
seurat_to_list <- function(obj){
 list(
   data = lapply(obj@assays,function(a){
     slots <- c("counts","data","scale.data")
     slots <- slots[sapply(slots,function(s){methods::.hasSlot(a,s)})]
     lapply(stats::setNames(slots,slots), function(s){
       methods::slot(a,s)
     })
   }) |> unlist(recursive = FALSE),
   obs = obj@meta.data,
   var = lapply(obj@assays,function(a){a@meta.features}),
   var_features = lapply(obj@assays,function(a){a@var.features}),
   reductions = lapply(obj@reductions,function(r){
     list(embedding = r@cell.embeddings,
          loadings = r@feature.loadings,
          loadings_projected = r@feature.loadings.projected)
   }),
   graphs = obj@graphs
 )
}
