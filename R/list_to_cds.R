#' Convert: \code{list} ==> \code{CellDataSet}
#'
#' @inheritParams converters
#' @inheritParams get_x
#' @returns A \code{CellDataSet} object.
#'
#' @export
#' @examples
#' obj <- example_obj("list")
#' obj2 <- list_to_cds(obj)
list_to_cds <- function(obj,
                        n = NULL,
                        merge_assays = TRUE,
                        verbose = TRUE,
                        ...){

  requireNamespace("monocle3")

  Xl <- get_x(obj = obj,
              n = n,
              simplify = FALSE,
              verbose = verbose,
              ...)
  var <- get_var(obj = obj,
                 n = n,
                 simplify = FALSE,
                 verbose = verbose,
                 ...)
  #### Check assays matche between Xl and var #####
  assay_nms <- get_assay_names(obj = Xl)
  if(!all(names(var) %in% assay_nms)){
    stopper("var must be a named list with all the same assays as the data:",
            paste(shQuote(assay_nms),collapse = " / "))
  }
  #### Create a set of CDS objects ####
  cds_list <- lapply(stats::setNames(names(Xl),
                                     names(Xl)),
                     function(nm){

    messager("Creating assay:",shQuote(nm),v=verbose)
    #### Get relevant var metadata ####
    var_i <- if(nm %in% names(var)){
      var[[nm]]
    } else if(names(var) %in% assay_nms){
      var[[names(var)[names(var) %in% assay_nms][1]]]
    } else {
      var[[1]]
    }
    #### Ensure same features ####
    var_i <- var_i[rownames(Xl[[nm]]),]
    #### Ensure var metadata ####
    if(!"gene_short_name" %in% names(var_i)){
      var_i$gene_short_name <- rownames(var_i)
    }
    #### Create CDS #####
    cds <- monocle3::new_cell_data_set(
      expression_data = Xl[[nm]],
      cell_metadata = get_obs(obj = obj,
                              verbose = verbose),
      gene_metadata = var_i
    )
    SummarizedExperiment::assayNames(cds) <- nm
    return(cds)
  })
  #### Find the most common assay dimensions ####
  n_rows <- as.integer(utils::tail(names(table(sapply(cds_list,nrow))),1))
  n_cols <- as.integer(utils::tail(names(table(sapply(cds_list,ncol))),1))
  for(nm in names(cds_list)){
    if((nrow(cds_list[[nm]])!=n_rows) ||
       (ncol(cds_list[[nm]])!=n_cols)){
      messager("Removing assay as its dimensions do not match the",
               "other assays:",shQuote(nm),v=verbose)
      cds_list[[nm]] <- NULL
    }
  }
  #### Merge into one object with multiple assays ####
  if(isTRUE(merge_assays)){

  } else {
    return(cds_list)
  }
}
