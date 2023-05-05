#' Convert: \code{list} ==> \code{AnnData}
#'
#' @export
#' @examples
#' obj <- example_obj("list")
#' obj2 <- list_to_anndata(obj)
list_to_anndata <- function(obj,
                            transpose=TRUE,
                            verbose=TRUE){
  # devoptera::args2vars(list_to_anndata)

  messager("+ list ==> AnnData",v=verbose)
  activate_conda(verbose=verbose)
  #### Get data ####
  ## Read in
  X <- get_data(obj = obj)
  ### Check if list
  if(is_list(X)){
    messager("Only first element in `data` will be used:",names(X)[1],
             v=verbose)
    X <- X[[1]]
  }
  X <- to_sparse(obj = X,
                 transpose = transpose)
  #### Get var ####
  ## Read in
  if(is.null(obj$var)){
    var <- data.frame(id=colnames(X),row.names = colnames(X))
  } else {
    var <- get_var(obj = obj,
                   verbose = verbose)
  }
  ## Check if list
  if(is_list(var)){
    messager("Only first element in `var` will be used:",names(var)[1],
             v=verbose)
    var <- var[[1]]
  }
  #### Get obs ####
  if(is.null(obj$obs)){
    obs <- data.frame(id=rownames(X), row.names = rownames(X))
  } else {
    obs <- get_obs(obj = obj,
                   verbose = verbose)
  }
  #### Get obsm ####
  obsm <- get_obsm(obj = obj,
                   verbose = verbose)
  #### Get varm ####
  varm <- get_varm(obj = obj,
                   verbose = verbose)
  varm <- check_varm(varm = varm,
                     X = X)
  #### Get uns ####
  uns <- get_uns(obj = obj,
                 verbose = verbose)
  #### Create new object ####
  messager("Constructing new AnnData object.",v=verbose)
  obj2 <- anndata::AnnData(
    X = X,
    obs = obs[rownames(X),,drop=FALSE],
    var = var[colnames(X),,drop=FALSE],
    obsm = obsm,
    varm = varm,
    uns = uns)
  return(obj2)
}
