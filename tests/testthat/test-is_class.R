test_that("is_class works", {

  dict <- class_dict()
  names(dict)
  #### ewce ####
  obj <- example_obj("ctd")
  testthat::expect_true(
    is_class(obj,"ewce")
  )
  #### matrix ####
  obj <- example_obj("matrix")
  testthat::expect_true(
    is_class(obj,"matrix")
  )
  #### loom ####
  library(Seurat)
  obj <- example_obj("loom")
  testthat::expect_true(
    is_class(obj,"loom")
  )
  #### se ####
  obj <- example_obj("se")
  testthat::expect_true(
    is_class(obj,"se")
  )
  obj <- example_obj("sce")
  testthat::expect_true(
    is_class(obj,"se")
  )
  obj <- example_obj("hdf5se")
  testthat::expect_true(
    is_class(obj,"se")
  )
  #### anndata ####
  obj <- example_obj("anndata")
  testthat::expect_true(
    is_class(obj,"anndata")
  )
  #### seurat ####
  obj <- example_obj("Seurat")
  testthat::expect_true(
    is_class(obj,"seurat")
  )
  #### h5Seurat ####
  obj <- example_obj("h5Seurat")
  testthat::expect_true(
    is_class(obj,"h5Seurat")
  )
  #### cds ####
  obj <- example_obj("cds")
  testthat::expect_true(
    is_class(obj,"cds")
  )
  #### cds ####
  obj <- example_obj("cds")
  testthat::expect_true(
    is_class(obj,"cds")
  )
})
