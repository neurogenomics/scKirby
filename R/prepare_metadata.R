#' Prepare metadata
#'
#' Standardise metadata for any single-cell object.
#' @export
#' @importFrom tidyr unite all_of
#' @examples
#' obj <- example_obj("anndata")
#' obj2 <- prepare_metadata()
prepare_metadata <- function(obj,
                             save_path = NULL,
                             species,
                             standardise_species = TRUE,
                             method = "gprofiler",
                             level_cols = list(),
                             extra_metadata = NULL,
                             return_metadata = FALSE,
                             verbose=TRUE){
  # devoptera::args2vars(prepare_metadata)

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
  if(isTRUE(return_metadata)){
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
