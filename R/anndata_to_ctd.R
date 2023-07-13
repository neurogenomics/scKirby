#' Convert: \code{AnnData} ==> \code{CellTypeDataset}
#'
#' @inheritParams orthogene::aggregate_rows
#'
#' @export
#' @import orthogene
#' @examples
#' obj <- example_obj("anndata")
#' obj2 <- anndata_to_ctd(obj, annotLevels=list(groups=NULL))
anndata_to_ctd <- function(obj,
                           annotLevels,
                           dataset = basename(tempfile()),
                           chunk_size = NULL,
                           agg_fun = "mean",
                           agg_method = c("monocle3",
                                          "stats"),
                           dropNA = TRUE,
                           standardise = TRUE,
                           as_sparse = TRUE,
                           as_DelayedArray = FALSE,
                           verbose = TRUE,
                           ...){
  # devoptera::args2vars(anndata_to_ctd)
  # annotLevels <- list(groups=NULL)

  #### Read data ####
  if(is.character(obj)){
    obj <- read_data(obj, backed="r")
  }
  if(isFALSE(obj$isbacked)){
    message("anndata object is not backed and may use more memory.")
  }
  #### Check annotLevels early ####
  annotLevels <- check_annotlevels(obj = obj,
                                   annotLevels = annotLevels)
  #### Aggregate cells over each level ####
  aggregate_rows <- utils::getFromNamespace("aggregate_rows","orthogene")
  rows <- nrow(obj)
  if(is.null(chunk_size)){
    chunk_size <- rows
  }
  chunks <- split(seq_len(rows), ceiling(seq_along(seq_len(rows))/chunk_size))
  ctd <- lapply(seq_len(length(annotLevels)), function(ix){
    messager("Generating CTD level:",ix,v=verbose)
    lvl <- annotLevels[[ix]]
    #### Aggregate within chunks ####
    X_list <- lapply(seq_len(length(chunks)), function(i) {
      messager("Processing chunk: ", i, "/",
               length(chunks),if (verbose > 1) "\n",
               parallel = TRUE,
               v = verbose)
      if (i == 1) verbose <- 2
      select <- as.integer(chunks[[i]])
      aggregate_rows(X = obj[select,]$X,
                                groupings = as.character(lvl[select]),
                                agg_fun = agg_fun,
                                agg_method = agg_method,
                                as_sparse = as_sparse,
                                as_DelayedArray = as_DelayedArray,
                                dropNA = dropNA,
                                verbose = verbose>1)
    })
    #### Aggregate across chunks ####
    X_agg <- aggregate_rows(
      X = do.call(what = rbind, X_list),
      groupings = unlist(lapply(X_list,rownames)),
      agg_fun = agg_fun,
      agg_method = agg_method,
      as_sparse = as_sparse,
      as_DelayedArray = as_DelayedArray,
      dropNA = dropNA,
      verbose = verbose) |> Matrix::t()
    remove(X_list)
    ctd_lvl <- EWCE::generate_celltype_data(
      exp = X_agg,
      groupName = paste0(dataset,"_level_",ix),
      annotLevels = list(colnames(X_agg)),
      return_ctd = TRUE,
      verbose = verbose)
    return(ctd_lvl$ctd[[1]])
  })
  #### Standardise CTD ####
  if(isTRUE(standardise)){
    ctd <- EWCE::standardise_ctd(ctd = ctd,
                                 dataset = dataset,
                                 verbose = verbose,
                                 ...)
  }
  return(ctd)
}
