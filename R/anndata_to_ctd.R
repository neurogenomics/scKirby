#' Convert: \code{AnnData} ==> \code{CellTypeDataset}
#'
#' @param standardise Run \link[EWCE]{standardise_ctd}.
#' @inheritParams converters
#' @inheritParams map_data
#' @inheritParams to_se
#' @inheritParams orthogene::aggregate_rows
#' @inheritParams EWCE::generate_celltype_data
#' @inheritParams EWCE::standardise_ctd
#' @inheritDotParams EWCE::standardise_ctd
#'
#' @export
#' @import orthogene
#' @importFrom utils getFromNamespace
#' @importFrom EWCE generate_celltype_data standardise_ctd
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
                           as_delayedarray = FALSE,
                           verbose = TRUE,
                           ...){
  # devoptera::args2vars(anndata_to_ctd)
  # annotLevels <- list(groups=NULL)

  messager_to_()
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
  chunks <- split(seq(rows), ceiling(seq_along(seq(rows))/chunk_size))
  ctd <- lapply(seq_len(length(annotLevels)), function(ix){
    messager("Generating CTD level:",ix,v=verbose)
    lvl <- annotLevels[[ix]]
    #### Aggregate within chunks ####
    X_list <- lapply(seq(length(chunks)), function(i) {
      messager("Processing chunk: ", i, "/",
               length(chunks),if (verbose > 1) "\n",
               parallel = TRUE,
               v = verbose)
      verbose <- if (i == 1) 2 else 1
      select <- as.integer(chunks[[i]])
      X <- aggregate_rows(X = obj[select,]$X,
                     groupings = as.character(lvl[select]),
                     agg_fun = agg_fun,
                     agg_method = agg_method,
                     as_sparse = as_sparse,
                     as_delayedarray = as_delayedarray,
                     dropNA = dropNA,
                     verbose = verbose>1)
      if(isTRUE(as_sparse)){
        X <- to_sparse(obj = X,
                       verbose = verbose)
      }
      return(X)
    })
    #### Aggregate across chunks ####
    X_agg <- aggregate_rows(
      X = do.call(what = rbind, X_list),
      groupings = unlist(lapply(X_list,rownames)),
      agg_fun = agg_fun,
      agg_method = agg_method,
      as_sparse = as_sparse,
      as_delayedarray = as_delayedarray,
      dropNA = dropNA,
      verbose = verbose) |> Matrix::t()
    ## Check for NAs
    X_agg <- fillna_sparse(X_agg)
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
