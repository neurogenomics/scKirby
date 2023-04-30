#' Convert: \code{list} ==> \code{AnnData}
#'
#' @export
#' @examples
#' obj <- example_obj("list")
#' obj2 <- list_to_anndata(obj)
list_to_anndata <- function(obj,
                            verbose=TRUE){
  # devoptera::args2vars(list_to_anndata)

  messager("+ list ==> AnnData",v=verbose)
  activate_conda(verbose=verbose)
  #### Get data ####
  if(length(obj$data)>1){
    messager("Only first assay being used:",names(obj$data)[1],v=verbose)
  }
  X <- to_sparse(obj =obj$data[[1]]) |>  Matrix::t()
  #### Get var ####
  if(is.null(obj$var)){
    var <- data.frame(id=colnames(X),row.names = colnames(X))
  } else if(!is.data.frame(obj$var) &&
            is.list(obj$var)){
    var <- obj$var[[1]]
  } else {
    var <- obj$var
  }
  #### Get obs ####
  if(is.null(obj$obs)){
    obs <- data.frame(id=rownames(X),row.names = rownames(X))
  } else {
    obs <- obj$obs
  }
  #### Get obsm ####
  obsm <- get_obsm(obj = obj,
                   verbose = verbose)
  #### Get varm ####
  varm <- get_varm(obj = obj,
                   verbose = verbose)
  varm <- check_varm(varm = varm,
                     X = X)
  #### Create new object ####
  messager("Constructing new AnnData object.",v=verbose)
  obj2 <- anndata::AnnData(
    X = X,
    obs = obs[rownames(X),,drop=FALSE],
    var = var[colnames(X),,drop=FALSE],
    obsm = obsm,
    varm = varm,
    uns = obj$uns)
  return(obj2)
}
