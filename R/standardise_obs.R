#' Prepare metadata
#'
#' Standardise metadata for any single-cell object.
#' @param standardise_species Standardise the name using
#' \link[orthogene]{map_species}.
#' @param level_cols Names of metadata columns that indicate the hierarchically
#' organized cell type annotation levels.
#' Each level can be a combination of other columns, e.g.:
#'  \code{
#'    list(level1="organ",
#'         level2=c("celltype"),
#'         level3=c("organ","celltype"))
#'    }
#' @param species The name of the species that the single-cell \code{obj}
#'  came from. Can be provided as the scientific name (e.g. "Mus musculus") or
#'  common name (e.g. "mouse").
#' @param save_path Path to save the modified \code{obj} to.
#' @param return_obs Return just standardised \code{obs}
#' instead of the entire \code{obj}.
#' @param extra_metadata Extra metadata to add to obs.
#'  Can be a:
#'  \itemize{
#'  \item{NULL : }{No extra metadata will be added (default).}
#'  \item{data.frame : }{With the same number of rows as observations (samples)
#'  in \code{obj}. NOTE: We assumed the rows of the \code{extra_metadata}
#'   data.frame are in the same order as the observations (samples)
#'  in \code{obj}}
#'  \item{list : }{A named list with a single entry per element,
#'  to be applied to all observation.}
#' }
#' @inheritParams converters
#' @inheritParams orthogene::map_species
#' @export
#' @importFrom tidyr unite all_of
#' @examples
#' obj <- example_obj()
#' obj2 <- standardise_obs()
standardise_obs <- function(obj,
                            save_path = NULL,
                            species,
                            standardise_species = TRUE,
                            method = "gprofiler",
                            level_cols = list(),
                            extra_metadata = NULL,
                            return_obs = FALSE,
                            verbose = TRUE){
  # devoptera::args2vars(standardise_obs)

  #### Extract metadata ####
  obs <- get_obs(obj = obj,
                 verbose = verbose)
  colnames(obs) <- tolower(colnames(obs))
  #### Add extra metadata ####
  for(col in extra_metadata){
    obs[[col]] <- extra_metadata[[col]]
  }
  #### Add levels ####
  for(lvl in names(level_cols)){
    messager("Adding column:",lvl,v=verbose)
    obs <- tidyr::unite(obs,
                        tidyr::all_of(tolower(level_cols[[lvl]])),
                        col = lvl,
                        sep = '.',
                        remove=FALSE)
    obs[lvl] <- obs$lvl
    obs$lvl <- NULL
  }
  #### Add species ####
  if(is.null(species)){
    species <- orthogene::infer_species(gene_df = obj$var_names,
                                        make_plot = FALSE,
                                        method = "homologene",
                                        verbose = verbose)$top_match
  }
  if(isTRUE(standardise_species)){
    species <- orthogene::map_species(species = species,
                                      method = method,
                                      verbose = verbose)
  }
  obs["species"] <- species
  #### Early return ####
  if(isTRUE(return_obs)){
    return(obs)
  }
  #### Add new obs ####
  obj <- set_obs(obj = obj,
                 obs = obs,
                 verbose = verbose)
  #### Make levels dictionary ####
  lvls_dict <- make_levels_dict(obs = obs,
                                level_cols = level_cols)
  obj <- set_uns(obj = obj,
                 uns = lvls_dict,
                 key = "levels_key",
                 verbose = verbose)
  #### Save ####
  if(!is.null(save_path)){
    save_path <- save_data(obj = obj,
                           save_path = save_path,
                           verbose = verbose)
  }
  return(obj)
}
