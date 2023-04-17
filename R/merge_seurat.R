#' Merge \link[Seurat]{Seurat} objects
#'
#' Merge a list of Seurat objects efficiently.
#' @param obj_list A named list of \link[Seurat]{Seurat} objects,
#' or paths to \link[Seurat]{Seurat} objects saved to disk as \emph{.rds} files.
#' @returns A merged \link[Seurat]{Seurat} object.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' obj_list <- Seurat::SplitObject(obj,split.by = "groups")
#' obj2 <- merge_seurat(obj_list)
merge_seurat <- function(obj_list,
                         batch_size = 5,
                         verbose = TRUE){

  #### Check names ####
  if(is.null(names(obj_list))){
    messager("Creating names for obj_list.",v=verbose)
    names(obj_list) <- paste0("obj",seq_len(length(obj_list)))
  }
  #### Create batches ####
  batches <- split(names(obj_list),
                   ceiling(seq_along( names(obj_list))/batch_size) )
  messager("Merging obj_list across",length(batches),"batch(es).",v=verbose)
  ##### Iterate over batches ####
  lapply(batches, function(b){
    lapply(b, function(nm){
      obj <- obj_list[[nm]]
      #### Import data #####
      if(is.character(obj) &&
         grepl("\\.rds$",obj,ignore.case = TRUE)){
        cat("Importing:",basename(obj),"\n")
        obj <- readRDS(obj)
      }
      #### Make cell names unique across objects ####
      obj@meta.data$cellid <- paste(nm,
                                    rownames(obj@meta.data),
                                    sep=".")
      Seurat::RenameCells(obj, new.names	= obj$cellid)
    }) |> Reduce(f=merge)
  }) |> Reduce(f=merge)
}
