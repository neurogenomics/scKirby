#' Convert matrix
#'
#' Convert a matrix object
#' across-species (gene orthologs) or within-species (gene synonyms).
#' @param obj A matrix object.
#' @inheritParams orthogene::aggregate_mapped_genes
#' @inheritParams orthogene::convert_orthologs
#' @inheritDotParams orthogene::convert_orthologs
#' @returns A nmatrix object.
#'
#' @keywords internal
#' @import orthogene
map_data_matrix <- function(obj,
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
                            as_delayedarray = FALSE,
                            sort_rows = TRUE,
                            test_species = NULL,
                            verbose = TRUE){
  assays <- map_data_assays(
    assays = list(obj=obj),
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
    as_delayedarray=as_delayedarray,
    sort_rows=sort_rows,
    test_species=test_species,
    verbose=verbose)
   return(assays[[1]])
}
