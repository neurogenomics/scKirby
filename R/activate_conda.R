#' Activate conda environment
#'
#' Activate a conda environment using a global variable.
#' Requires the \code{echoconda} package (optional dependency).
#' If \code{echoconda} is not installed, this function does nothing.
#' @param conda_env Name of the conda environment to activate.
#' @param method Method for activating the environment (passed to echoconda).
#' @param verbose Print messages.
#' @returns Null (invisibly).
#'
#' @export
#' @examples
#' \dontrun{
#' activate_conda()
#' }
activate_conda <- function(conda_env = Sys.getenv("CONDA_ENV"),
                           method = "reticulate",
                           verbose = TRUE){

  #### Check if echoconda is available ####
  if(!requireNamespace("echoconda", quietly = TRUE)){
    messager("echoconda not installed. Skipping conda activation.",
             "Install with: remotes::install_github('RajLabMSSM/echoconda')",
             v = verbose)
    return(invisible(NULL))
  }
  #### Set default conda env ####
  if(conda_env=="") conda_env <- "r-reticulate"
  #### Only activate if not already activated ####
  ## Use get() to avoid parse-time namespace lookup
  which_env <- get("which_env", asNamespace("echoconda"))
  activate_env <- get("activate_env", asNamespace("echoconda"))
  if(which_env()!=conda_env){
    activate_env(conda_env = conda_env,
                 method = method,
                 verbose = verbose)
  }
  return(invisible(NULL))
}
