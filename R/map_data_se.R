#' Convert \link[SummarizedExperiment]{SummarizedExperiment}
#'
#' Convert a \link[SummarizedExperiment]{SummarizedExperiment} object
#' across-species (gene orthologs) or within-species (gene synonyms).
#' Also works for \code{SingleCellExperiment} objects, or any other extensions
#' of the \link[SummarizedExperiment]{SummarizedExperiment} class.
#' @param obj A \link[SummarizedExperiment]{SummarizedExperiment} object.
#' @inheritParams orthogene::aggregate_mapped_genes
#' @inheritParams orthogene::convert_orthologs
#' @inheritDotParams orthogene::convert_orthologs
#' @returns \link[SummarizedExperiment]{SummarizedExperiment}
#'
#' @keywords internal
#' @import orthogene
map_data_se <- function(obj,
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
                        as_sparse = TRUE,
                        as_DelayedArray = FALSE,
                        sort_rows = TRUE,
                        test_species = NULL,
                        verbose = TRUE){
    # devoptera::args2vars(map_data_se)
    # obj <- example_obj("se")

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
      verbose=verbose)
    #### Construct row data using gene map ####
    rd <- map_data_rowdata(
      genes = rownames(assays[[1]]),
      original_rowdata = SummarizedExperiment::rowData(obj))
    #### Construct new SummarizedExperiment ####
    obj2 <- SummarizedExperiment::SummarizedExperiment(
        assays = assays,
        rowData = rd[rownames(assays[[1]]),,drop=FALSE],
        colData = obj@colData,
        metadata = obj@metadata)
    #### Convert back to originla subclass ###
    if(methods::is(obj,"SingleCellExperiment")){
      obj2 <- se_to_sce(obj = obj2, verbose = FALSE)
    }
    return(obj2)
}
