#' Convert \link[SeuratObject]{Seurat}
#'
#' Convert a \link[SeuratObject]{Seurat} object
#' across-species (gene orthologs) or within-species (gene synonyms).
#' @param obj A \link[SeuratObject]{Seurat} object.
#' @inheritParams orthogene::aggregate_mapped_genes
#' @inheritParams orthogene::convert_orthologs
#' @inheritDotParams orthogene::convert_orthologs
#' @returns \link[SeuratObject]{Seurat}
#'
#' @keywords internal
#' @import orthogene
map_data_seurat <- function(obj,
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
                            sort_rows = TRUE,
                            test_species = NULL,
                            verbose = TRUE){
  # devoptera::args2vars(map_data_seurat)
  # obj <- example_obj("seurat")

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
    as_DelayedArray=as_DelayedArray,
    sort_rows=sort_rows,
    test_species=test_species,
    verbose=verbose)
  #### Construct list of Assay objects #####
  aobj_list <- lapply(stats::setNames(Seurat::Assays(obj),
                                      Seurat::Assays(obj)),
                      function(nm){
    prefix <- paste0("^",nm,"\\.")
    a1 <- assays[grep(prefix, names(assays))]
    names(a1) <- gsub(prefix,"",names(a1))
    aobj <- SeuratObject::CreateAssayObject(counts = a1$counts,)
    if("scale.data" %in% names(a1)){
      aobj@scale.data <- a1$scale.data
    }
    if("data" %in% names(a1)){
      aobj@data <- a1$data
    }
    #### Construct row data using gene map ####
    rd <- map_data_rowdata(
      genes = rownames(aobj@counts),
      original_rowdata = obj@assays[[nm]]@meta.features)
    aobj@meta.features <- rd
    return(aobj)
  })
  #### Construct new SeuratObject ####
  obj2 <- SeuratObject::CreateSeuratObject(
    counts = aobj_list[[1]],
    assay = names(aobj_list)[[1]],
    meta.data = obj@meta.data)
  #### Add extra assays ####
  if(length(aobj_list)>1){
    for(nm in names(aobj_list)[-1])
    obj2[[nm]] <- aobj_list[[nm]]
  }
  return(obj2)
}
