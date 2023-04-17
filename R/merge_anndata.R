#' Merge \link[anndata]{AnnData} objects
#'
#' Merge a list of Seurat objects efficiently.
#' @param obj_list A named list of \link[anndata]{AnnData} objects,
#' or paths to \link[anndata]{AnnData} objects
#' saved to disk as \emph{.h5ad} files.
#' @inheritParams anndata::concat
#' @returns A merged \link[anndata]{AnnData} object.
#'
#' @export
#' @examples
#' obj <- example_obj("anndata")
#' obj_list <- Seurat::SplitObject(obj,split.by = "groups")
#' obj2 <- merge_seurat(obj_list)
merge_anndata <- function(obj_list,
                          axis = 0L,
                          join = "outer",
                          merge = NULL,
                          uns_merge = NULL,
                          label = NULL,
                          keys = NULL,
                          index_unique = NULL,
                          fill_value = NULL,
                          pairwise = FALSE,
                          verbose = TRUE){

  #### Check names ####
  if(is.null(names(obj_list))){
    messager("Creating names for obj_list.",v=verbose)
    names(obj_list) <- paste0("obj",seq_len(length(obj_list)))
  }
  #### Create batches ####
  messager("Merging obj_list.",v=verbose)
  ##### Iterate objs ####
  lapply(names(obj_list), function(nm){
    obj <- obj_list[[nm]]
    #### Import data #####
    if(is.character(obj) &&
       grepl("\\.h5ad$",obj,ignore.case = TRUE)){
      cat("Importing:",basename(obj),"\n")
      obj <- anndata::read_h5ad(obj)
    }
    #### Make cell names unique across objects ####
    obj$obs$cellid <- paste(nm,
                            rownames(obj$obs),
                            sep=".")
    obj$obs_names <- obj$obs$cellid
    obj$obs_names_make_unique()
    return(obj)
  }) |> anndata::concat(axis = axis,
                        join = join,
                        merge = merge,
                        uns_merge = uns_merge,
                        label = label,
                        keys = keys,
                        index_unique = index_unique,
                        fill_value = fill_value,
                        pairwise = pairwise)
}
