#' Compute all pairwise trait correlations
#'
#' Computes pairwise correlations for all traits
#'  using either \link[WGCNA]{cor} (faster) or
#' \link[stats]{cor} (slower).
#' @inheritParams get_obsm
#' @inheritParams stats::cor
#' @inheritParams WGCNA::cor
#' @returns A single-cell object with a graph.
#'
#' @keywords internal
calc_cor <- function(obj,
                     keys = NULL,
                     graph_key = NULL,
                     assay = NULL,
                     layer = NULL,
                     transpose = FALSE,
                     method = "pearson",
                     fill_na = NULL,
                     use = "all.obs",
                     nThreads = 1,
                     verbose = TRUE) {
    # devoptera::args2vars(calc_cor)

  #### Extract relevant matrix ####
  if (!is.null(keys)) {
      X <- get_obsm(
          obj = obj,
          keys = keys,
          n = 1,
          verbose = verbose)
      #### Important! must transpose ####
      X <- to_sparse(obj = X, transpose = TRUE)
  } else {
      X <- get_x(
          obj = obj,
          assay = assay,
          layer = layer,
          transpose = transpose,
          n = 1,
          verbose = verbose)
  }
  if (!is.null(fill_na)) X[is.na(X)] <- fill_na
  #### Compute corr ####
  if (is_installed(pkg = "WGCNA")) {
      messager("Computing r with WGCNA.", v = verbose)
      Xcor <- WGCNA::cor(x = X,
                         method = method,
                         use = use,
                         nThreads = nThreads,
                         verbose = verbose)
  } else {
      messager("Computing r with stats.", v = verbose)
      Xcor <- stats::cor(x = X,
                         use = use,
                         method = method)
  }
  return(Xcor)
}
