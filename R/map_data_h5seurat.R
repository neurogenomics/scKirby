#' Convert \code{h5Seurat}
#'
#' Convert a \link[SeuratDisk]{h5Seurat} object
#' across-species (gene orthologs) or within-species (gene synonyms).
#' @param obj A \link[SeuratDisk]{h5Seurat} object.
#' @inheritParams orthogene::aggregate_mapped_genes
#' @inheritParams orthogene::convert_orthologs
#' @inheritDotParams orthogene::convert_orthologs
#' @returns A \link[SeuratDisk]{h5Seurat} object.
#'
#' @keywords internal
#' @import orthogene
map_data_h5seurat <- function(obj,
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
                            as_delayedarray = FALSE,
                            sort_rows = TRUE,
                            test_species = NULL,
                            verbose = TRUE){
  # devoptera::args2vars(map_data_h5seurat)
  # obj <- example_obj("h5seurat")


  # obj[["assays/RNA"]]$ls()
  # hd_data <- obj[["assays/RNA/data"]]
  # SeuratDisk::as.list(hd_data)
  # X <- SeuratObject::as.sparse(hd_data)[1:100,1:100, drop=FALSE]
  # var <- as.data.frame(obj[["assays/RNA/meta.features"]])
  # rownames(X) <- rownames(var)
  # SeuratDisk::as.list(obj[["assays"]])
  # obj

# SeuratObject::as.Seurat(, )


  l <- seurat_to_list(obj = obj)
  assays <- map_data_assays(
    assays = l$data,
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
  #### Construct new SeuratObject ####
  l$data <- assays
  obj2 <- list_to_seurat(obj = l,
                         verbose = verbose)
  return(obj2)
}
