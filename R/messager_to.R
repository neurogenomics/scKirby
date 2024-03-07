#' Messager to
#'
#' Automated messager that automatically infers the function from which it's
#'  being called.
#' @keywords internal
messager_to <- function(element=c("name","link")){

  envir <- parent.frame(n = 2L)
  v <- if("verbose" %in% ls(envir = envir)){
    get("verbose",envir = envir)
  } else {
    TRUE
  }
  element <- element[1]
  func <- names(parent.frame(2L))
  splt <- strsplit(func,"_")[[1]]

  d <- dict_acronym()
  messager("Converting:",
           d[[splt[1]]][[element]],
           "==>",
           d[[splt[3]]][[element]],
           v=v)
}
