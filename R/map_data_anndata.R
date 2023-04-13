#' Convert \link[anndata]{AnnData}
#'
#' Convert an \link[anndata]{AnnData} object
#' across-species (gene orthologs) or within-species (gene synonyms).
#' @inheritParams orthogene::aggregate_mapped_genes
#' @inheritParams orthogene::convert_orthologs
#' @inheritParams echoconda::activate_env
#' @returns \link[anndata]{AnnData}
#'
#' @keywords internal
#' @importFrom echoconda activate_env
map_data_anndata <- function(obj,
                             gene_map = NULL,
                             input_col = "input_gene",
                             output_col = "ortholog_gene",
                             standardise_genes = FALSE,
                             input_species = NULL,
                             output_species = input_species,
                             method = c(
                               "homologene",
                               "gprofiler",
                               "babelgene"
                             ),
                             drop_nonorths = TRUE,
                             non121_strategy =
                               "drop_both_species",
                             agg_fun = NULL,
                             mthreshold = Inf,
                             as_sparse = FALSE,
                             as_DelayedArray = FALSE,
                             sort_rows = FALSE,
                             test_species = NULL,
                             conda_env = "r-reticulate",
                             verbose = TRUE){
  # devoptera::args2vars(map_data_anndata)
  # obj <- example_obj("anndata")

  #### Activate conda env with anndata installed ####
  echoconda::activate_env(conda_env = conda_env,
                          method = "reticulate",
                          verbose = verbose)

  assays <- list(X=obj$X,
                 raw=obj$raw)
  #### Convert and transpose ####
  assays <- lapply(assays, function(a){
    if(!is.null(a)){
      if(methods::is(a,"RawR6")){
        a <- as.matrix(a)
      }
      Matrix::t(a)
    }
  })
  assays <- map_data_assays(
    assays = assays,
    gene_map=gene_map,
    input_col=input_col,
    output_col=output_col,
    standardise_genes=standardise_genes,
    input_species=input_species,
    output_species=output_species,
    method=method,
    drop_nonorths=drop_nonorths,
    non121_strategy=non121_strategy,
    agg_fun=agg_fun,
    mthreshold=mthreshold,
    as_sparse=as_sparse,
    as_DelayedArray=as_DelayedArray,
    sort_rows=sort_rows,
    test_species=test_species,
    verbose=verbose)
  #### Construct row data using gene map ####
  rd <- map_data_rowdata(
    genes = rownames(assays[[1]]),
    original_rowdata = obj$var)
  #### Construct new SummarizedExperiment ####
  obj2 <- anndata::AnnData(
    X = Matrix::t(assays$X),
    raw = if(!is.null(assays$raw))Matrix::t(assays$raw),
    obs = obj$obs,
    obsm = obj$obsm,
    obsp = obj$obsp,
    var = rd,
    varm = obj$varm,
    varp = obj$varp,
    uns = obj$uns,
    # layers = obj$layers, ## OMIT!: Causes "Error: KeyError: 1"
    filename = obj$filename)
  return(obj2)
}
