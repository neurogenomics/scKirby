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
  #### seurat ####
  obj <- example_obj("seurat")
  obj2 <- map_data(obj = obj,
                   input_species = "human",
                   output_species = "mouse")
  testthat::expect_lt(nrow(obj2),
                      nrow(obj))
})
