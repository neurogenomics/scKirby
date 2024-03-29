#' Merge \code{Seurat} objects
#'
#' Merge a list of Seurat objects efficiently.
#' @param obj_list A named list of \link[SeuratObject]{Seurat} objects,
#' or paths to \link[SeuratObject]{Seurat} objects saved to disk
#' as \emph{.rds} files.
#' @inheritParams map_data
#' @returns A merged \link[SeuratObject]{Seurat} object.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj_list <- Seurat::SplitObject(obj,split.by = "groups")
#' obj2 <- merge_seurat(obj_list)
merge_seurat <- function(obj_list,
                         chunk_size = 5,
                         verbose = TRUE){

  #### Check names ####
  if(is.null(names(obj_list))){
    messager("Creating names for obj_list.",v=verbose)
    names(obj_list) <- paste0("obj",seq_len(length(obj_list)))
  }
  #### Create batches ####
  batches <- split(names(obj_list),
                   ceiling(seq_along( names(obj_list))/chunk_size) )
  messager("Merging obj_list across",length(batches),"batch(es).",v=verbose)
  ##### Iterate over batches ####
  lapply(batches, function(b){
    lapply(b, function(nm){
      obj <- obj_list[[nm]]
      #### Import data #####
      if(is.character(obj)){
        obj <- read_data(path = obj,
                         verbose = verbose)
      }
      #### Make cell names unique across objects ####
      obj@meta.data$cellid <- paste(nm,
                                    rownames(obj@meta.data),
                                    sep=".")
      Seurat::RenameCells(obj, new.names	= obj$cellid)
    }) |> Reduce(f=merge)
  }) |> Reduce(f=merge)
}
