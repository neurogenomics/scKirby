test_that("to_seurat works from list", {
  obj <- example_obj("list")
  obj2 <- to_seurat(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "seurat"))
})

test_that("to_seurat works from matrix", {
  obj <- example_obj("matrix")
  obj2 <- to_seurat(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "seurat"))
})

test_that("to_seurat works from SE", {
  obj <- example_obj("se")
  obj2 <- to_seurat(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "seurat"))
})

test_that("to_seurat works from SCE", {
  obj <- example_obj("sce")
  obj2 <- to_seurat(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "seurat"))
})

test_that("to_seurat returns seurat object as-is", {
  obj <- example_obj("seurat")
  obj2 <- to_seurat(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "seurat"))
})

test_that("to_seurat works from CTD", {
  obj <- example_obj("ctd")
  obj2 <- to_seurat(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "seurat"))
})

test_that("to_seurat works from h5Seurat", {
  skip_if_no_seuratdisk()
  obj <- example_obj("h5seurat")
  obj2 <- to_seurat(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "seurat"))
})

test_that("to_seurat works from loom", {
  skip_if_no_loomr()
  library(Seurat)
  obj <- example_obj("loom")
  obj2 <- to_seurat(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "seurat"))
})

test_that("to_seurat works from anndata", {
  skip_if_no_anndata()
  obj <- example_obj("anndata")
  obj2 <- to_seurat(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "seurat"))
})

test_that("to_seurat works from CDS", {
  skip_if_no_monocle3()
  obj <- example_obj("cds")
  obj2 <- to_seurat(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "seurat"))
})
