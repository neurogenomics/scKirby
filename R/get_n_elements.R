#' Get N elements
#'
#' Get only the first N elements from a named list.
#' @param l A list.
#' @param n The maximum number of elements to return from that list.
#' @param simplify When the final length of the list is 1, unlist the object.
#' @param verbose Print messages.
#' @returns A shortened list or a single element.
#'
#' @export
#' @examples
#' l <- lapply(seq(4),prod,3)
#' l1 <- get_n_elements(l, n=1)
#' l2 <- get_n_elements(l, n=2, simplify = FALSE)
get_n_elements <- function(l,
                           n=1,
                           simplify=TRUE,
                           verbose=TRUE){

  #### Check if list ####
  vname <-  deparse(substitute(l))
  if(!is.list(l) ||
     is_class(l,"matrix")){
    # messager(vname,"is not a list. Returning unmodified object.",v=verbose)
    return(l)
  }
  #### Check number of elements ####
  if(is.null(n)) return(l)
  #### Get names ####
  if(is.null(names(l))){
    names(l) <- as.character(seq(length(l)))
  }
  keys <- names(l)
  #### Check number of elements requested ####
  if(length(n)>length(l)){
    messager("n cannot be greater than the length of l.",
             paste0("Setting n=",n),v=verbose)
    n <- length(l)
  }
  #### £xtract elements ####
  if(is.list(l) &&
     length(l)>n){
      messager(paste0(">",n),shQuote(vname),"element identified.",
               "Using first",n,"element(s) only:",
               paste0(if(n>1) "\n - " else '',
                     shQuote(keys[seq(n)]),collapse = ""),
               v=verbose)

  }
  l <- l[seq(n)]
  #### Simplify ####
  if(isTRUE(simplify)){
    if(length(l)==1){
      messager("Simplifying list",paste0(vname,"."),v=verbose)
      l <- l[[1]]
    } else {
      messager("List cannot be simplified when n>1.",v=verbose)
    }
  }
  #### Return ####
  return(l)
}
