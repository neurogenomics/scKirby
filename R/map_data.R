#' Convert single-cell data
#'
#' Convert a single-cell data object
#' across-species (gene orthologs) or within-species (gene synonyms).
#' @param obj A single-cell data object belonging to one of the
#' following classes:
#' \itemize{
#' \item{\link[SummarizedExperiment]{SummarizedExperiment}}
#' \item{\link[SingleCellExperiment]{SingleCellExperiment}}
#' \item{\link[SeuratObject]{Seurat}}
#' \item{\link[AnnData]{AnnData}}
#' \item{
#' \link[Matrix]{Matrix} or
#' \link[data.frame]{data.frame} or
#' \link[data.table]{data.table}
#' }
#' }
#' @param chunk_size An integer indicating number of cells to include per chunk.
#' This can be a more memory efficient and scalable way of aggregating on-disk
#' data formats like \link[AnnData]{AnnData},
#' rather than reading in the entire matrix into memory at once
#'  (default: \code{NULL}).
#' @inheritParams orthogene::infer_species
#' @inheritParams orthogene::convert_orthologs
#' @inheritParams orthogene::aggregate_mapped_genes
#' @returns An aggregated/expanded version of the input single-cell data object.
#'
#' @export
#' @import orthogene
#' @examples
#' obj <- example_obj("ad")
#' obj2 <- map_data(obj = obj,
#'                  input_species = "human",
#'                  output_species = "mouse")
map_data <- function(obj,
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
                     non121_strategy = "drop_both_species",
                     agg_fun = "sum",
                     mthreshold = Inf,
                     as_sparse = FALSE,
                     as_DelayedArray = FALSE,
                     sort_rows = FALSE,
                     test_species = NULL,
                     chunk_size = NULL,
                     verbose = TRUE,
                     ...){
  # devoptera::args2vars(map_data)
  # obj <- TabulaMurisData::TabulaMurisSmartSeq2()
  # obj <- anndata::read_h5ad("~/Downloads/EarthwormCellLandscape.h5ad")
  # gene_map <- orthogene:::earthworm2human_map()
  # input_col <- "qseqid"; output_col <- "sseqid"

  #### Infer species ####
  if(is.null(input_species)){
      input_species <- orthogene::infer_species(
        gene_df = rownames(obj),
        make_plot = FALSE,
        method = method,
        test_species = test_species,
        verbose = verbose)$top_match
      if(is.null(output_species)) output_species <- input_species
  }
  if(input_species==output_species){
    messager("Setting non121_strategy='keep_both_species'",
             "when input_species==output_species.",v=verbose)
    non121_strategy <- "keep_both_species"
  }
  if(is_class(obj,"se")){
      obj2 <- map_data_se(
          obj = obj,
          gene_map = gene_map,
          input_col = input_col,
          output_col = output_col,
          standardise_genes = standardise_genes,
          input_species = input_species,
          output_species = output_species,
          method = method,
          drop_nonorths = drop_nonorths,
          non121_strategy = non121_strategy,
          agg_fun = agg_fun,
          mthreshold = mthreshold,
          as_sparse = as_sparse,
          as_DelayedArray = as_DelayedArray,
          sort_rows = sort_rows,
          test_species = test_species,
          verbose = verbose)
  } else if(is_class(obj,"anndata")){
      obj2 <- map_data_anndata(
          obj = obj,
          gene_map = gene_map,
          input_col = input_col,
          output_col = output_col,
          standardise_genes = standardise_genes,
          input_species = input_species,
          output_species = output_species,
          method = method,
          drop_nonorths = drop_nonorths,
          non121_strategy = non121_strategy,
          agg_fun = agg_fun,
          mthreshold = mthreshold,
          as_sparse = as_sparse,
          as_DelayedArray = as_DelayedArray,
          sort_rows = sort_rows,
          test_species = test_species,
          chunk_size = chunk_size,
          verbose = verbose)
  } else if(is_class(obj,"seurat")){
    obj2 <- map_data_seurat(
      obj = obj,
      gene_map = gene_map,
      input_col = input_col,
      output_col = output_col,
      standardise_genes = standardise_genes,
      input_species = input_species,
      output_species = output_species,
      method = method,
      drop_nonorths = drop_nonorths,
      non121_strategy = non121_strategy,
      agg_fun = agg_fun,
      mthreshold = mthreshold,
      as_sparse = as_sparse,
      as_DelayedArray = as_DelayedArray,
      sort_rows = sort_rows,
      test_species = test_species,
      verbose = verbose)
  } else if(is_class(obj,"h5seurat")){
    obj2 <- map_data_h5seurat(
      obj = obj,
      gene_map = gene_map,
      input_col = input_col,
      output_col = output_col,
      standardise_genes = standardise_genes,
      input_species = input_species,
      output_species = output_species,
      method = method,
      drop_nonorths = drop_nonorths,
      non121_strategy = non121_strategy,
      agg_fun = agg_fun,
      mthreshold = mthreshold,
      as_sparse = as_sparse,
      as_DelayedArray = as_DelayedArray,
      sort_rows = sort_rows,
      test_species = test_species,
      verbose = verbose)
  } else {
    stopper("Unsupported obj class.")
  }
  return(obj2)
}
