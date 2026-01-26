test_that("is_class works for basic types", {
  dict <- dict_class()
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
  #### seurat ####
  obj <- example_obj("Seurat")
  testthat::expect_true(
    is_class(obj,"seurat")
  )
})

test_that("is_class works for loom", {
  skip_if_no_loomr()
  library(Seurat)
  obj <- example_obj("loom")
  testthat::expect_true(
    is_class(obj,"loom")
  )
})

test_that("is_class works for anndata", {
  skip_if_no_anndata()
  obj <- example_obj("anndata")
  testthat::expect_true(
    is_class(obj,"anndata")
  )
})

test_that("is_class works for h5Seurat", {
  skip_if_no_seuratdisk()
  obj <- example_obj("h5Seurat")
  testthat::expect_true(
    is_class(obj,"h5Seurat")
  )
})

test_that("is_class works for cds", {
  skip_if_no_monocle3()
  obj <- example_obj("cds")
  testthat::expect_true(
    is_class(obj,"cds")
  )
})
