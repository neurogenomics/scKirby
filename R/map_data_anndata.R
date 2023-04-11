#' Convert \link[anndata]{AnnData}
#'
#' Convert an \link[anndata]{AnnData} object
#' across-species (gene orthologs) or within-species (gene synonyms).
#' @inheritParams orthogene::aggregate_mapped_genes
#' @inheritParams orthogene::convert_orthologs
#' @returns \link[anndata]{AnnData}
#'
#' @keywords internal
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
                             verbose = TRUE,
                             ...){
  # devoptera::args2vars(map_data_anndata)
  # obj <- example_obj("anndata")

  assays <- map_data_assays(
    assays = SummarizedExperiment::assays(obj),
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
    verbose=verbose,
    ...)
  #### Construct row data using gene map ####
  genes <- rownames(assays[[1]])
  rd <- data.frame(input_gene=names(genes),
                   ortholog_gene=unname(genes)) |>
    #### merge with the original metadata ####
  merge(y = SummarizedExperiment::rowData(obj),
        all.x = TRUE,
        by.x = "input_gene",
        by.y = 0)
  rownames(rd) <- rd$ortholog_gene
  #### Construct new SummarizedExperiment ####
  obj2 <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    rowData = rd[rownames(assays[[1]]),],
    colData = obj@colData,
    metadata = obj@metadata)
  return(obj2)
}
