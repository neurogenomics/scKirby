#' Activate conda environment
#'
#' Activate a conda environment using a global variable.
#' @returns Null
#' @inheritParams echoconda::activate_env
#'
#' @export
#' @examples
#' activate_conda()
activate_conda <- function(conda_env = Sys.getenv("CONDA_ENV"),
                           method = "reticulate",
                           verbose = TRUE){

  #### Set default conda env ####
  if(conda_env=="") conda_env <- "r-reticulate"
  #### Only activate if not already activated ####
  if(echoconda::which_env()!=conda_env){
    echoconda::activate_env(conda_env = conda_env,
                            method = method,
                            verbose = verbose)
  }
}
