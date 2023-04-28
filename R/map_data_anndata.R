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
#' @importFrom orthogene map_orthologs
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
                             agg_fun = "sum",
                             mthreshold = Inf,
                             as_sparse = FALSE,
                             as_DelayedArray = FALSE,
                             sort_rows = FALSE,
                             test_species = NULL,
                             conda_env = "r-reticulate",
                             chunk_size = NULL,
                             verbose = TRUE){
  # devoptera::args2vars(map_data_anndata)
  # obj <- example_obj("anndata")

  #### Activate conda env with anndata installed ####
  echoconda::activate_env(conda_env = conda_env,
                          method = "reticulate",
                          verbose = verbose)
  #### Convert orthologs ####
  ## Use the same map for each chunk to reduce API queries.
  gene_map  <- orthogene::map_orthologs(genes = obj$var_names,
                                        gene_map = gene_map,
                                        method = method,
                                        input_species = input_species,
                                        output_species = output_species,
                                        standardise_genes = standardise_genes,
                                        mthreshold = mthreshold,
                                        input_col = input_col,
                                        output_col = output_col,
                                        verbose = verbose)
  #### Split matrix into chunks and process each one ####
  ## Improves memory efficiency
  assays <- map_data_anndata_chunked(obj=obj,
                                     chunk_size=chunk_size,
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
  messager("Reconstructing anndata object.",v=verbose)
  #### Construct row data using gene map ####
  rd <- map_data_rowdata(
    genes = rownames(assays[[1]]),
    original_rowdata = obj$var)
  #### Construct new SummarizedExperiment ####
  ## Remove PC slot to avoid error:
  ## Error: ValueError: Value passed for key 'PCs' is of incorrect shape.
  obj$obsm$X_pca <- NULL
  obj$varm$PCs <- NULL
  obj2 <-  anndata::AnnData(
    X = Matrix::t(assays$X),
    ### OMIT! 'raw': Causes Error: ValueError: The truth value of an array with more
    ## than one element is ambiguous. Use a.any() or a.all().
    # raw = if(!is.null(assays$raw))Matrix::t(assays$raw),
    obs = obj$obs,
    var = rd[rownames(assays$X),],
    obsm = obj$obsm,
    obsp = obj$obsp,
    ### OMIT! 'varm': Causes error due to mismatches between old and new vars
    # varm = obj$varm,
    varp = obj$varp,
    uns = obj$uns
    ## OMIT! 'layers': Causes "Error: KeyError: 1"
    # layers = obj$layers,
    ## OMIT! 'filename' can cause issues
    # filename = obj$filename
    )
  return(obj2)
}
