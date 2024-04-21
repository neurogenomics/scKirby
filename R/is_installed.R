#' Is package installed
#'
#' @param pkg Name of an R package to check.
#' @param error If \code{TRUE}, produces an error if the \code{pkg}
#' isn't installed.
#' If \code{FALSE}, produces a warning if the package isn't installed.
#' @returns Boolean
#'
#' @keywords internal
is_installed <- function(pkg,
                         error=FALSE) {
  bool <- requireNamespace(pkg, quietly = TRUE)
  func <- if(error) stop else warning
  if(!bool){
    txt <- paste("Warning: Must install",shQuote(pkg),"to use this feature.")
    func(txt)
  }
  return(bool)
}
