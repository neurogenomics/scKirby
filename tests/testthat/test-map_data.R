test_that("map_data works", {

  #### SummarizedExperiment ####
  obj <- example_obj("se")
  obj2 <- map_data(obj = obj,
                   input_species = "human",
                   output_species = "mouse")
  testthat::expect_lt(nrow(obj2),
                      nrow(obj))
  #### SingleCellExperiment ####
  obj <- example_obj("sce")
  obj2 <- map_data(obj = obj,
                   input_species = "human",
                   output_species = "mouse")
  testthat::expect_lt(nrow(obj2),
                      nrow(obj))
  #### anndata ####
  obj <- example_obj("ad")
  obj2 <- map_data(obj = obj,
                   input_species = "human",
                   output_species = "mouse")
  testthat::expect_lt(ncol(obj2),
                      ncol(obj))
  testthat::expect_false(is(obj2$X,"dgCMatrix"))
  #### anndata: chunked ####
  obj <- example_obj("ad")
  obj2 <- map_data(obj = obj,
                   input_species = "human",
                   output_species = "mouse",
                   chunk_size = 10,
                   as_sparse = TRUE)
  testthat::expect_lt(ncol(obj2),
                      ncol(obj))
  testthat::expect_true(is(obj2$X,"dgCMatrix"))
  #### seurat ####
  obj <- example_obj("seurat")
  obj2 <- map_data(obj = obj,
                   input_species = "human",
                   output_species = "mouse")
  testthat::expect_lt(nrow(obj2),
                      nrow(obj))


  #### anndata: using gene_map across species ####
  obj <- example_obj("ad")
  gene_map <- orthogene::map_orthologs(genes = colnames(obj),
                                       input_species = "human",
                                       output_species = "mouse")
  obj2 <- map_data(obj = obj,
                   gene_map = gene_map,
                   input_species = "human",
                   output_species = "mouse")
  testthat::expect_lt(ncol(obj2),
                      ncol(obj))

  #### anndata: using gene_map within species ####
  ### Simply renames the genes if the number is equal
  obj <- example_obj("ad")
  gene_map <- data.frame(input=paste("gene",
                                     sort(rep(seq_len(ncol(obj)/2),2)),c(1,2),
                                     sep="."),
                         output=colnames(obj))
  obj$var_names <- gene_map$input
  obj2 <- map_data(obj = obj,
                   gene_map = gene_map,
                   input_species = "human",
                   output_species = "human",
                   input_col = "input",
                   output_col = "output")
  testthat::expect_equal(ncol(obj2),
                         ncol(obj))
  ### Aggregate if not equal
  obj <- example_obj("ad")
  gene_map <- data.frame(input=colnames(obj),
                         output=sort(rep(paste0("gene",seq_len(ncol(obj)/2)),2)))
  obj2 <- map_data(obj = obj,
                   gene_map = gene_map,
                   input_species = "human",
                   output_species = "human",
                   input_col = "input",
                   output_col = "output")
  testthat::expect_lte(ncol(obj2),
                       ncol(obj))

})
