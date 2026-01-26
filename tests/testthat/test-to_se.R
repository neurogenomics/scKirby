test_that("to_se works from Seurat", {
  obj <- example_obj("seurat")
  obj2 <- to_se(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "se"))
})

test_that("to_se works with as_sce=TRUE", {
  obj <- example_obj("seurat")
  obj2 <- to_se(obj, as_sce = TRUE, verbose = FALSE)
  testthat::expect_true(methods::is(obj2, "SingleCellExperiment"))
})

test_that("to_se works from matrix", {
  obj <- example_obj("matrix")
  obj2 <- to_se(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "se"))
})

test_that("to_se works from CTD", {
  obj <- example_obj("ctd")
  obj2 <- to_se(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "se"))
})

test_that("to_se returns SE object as-is", {
  obj <- example_obj("se")
  obj2 <- to_se(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "se"))
})

test_that("to_se works from anndata", {
  skip_if_no_anndata()
  obj <- example_obj("anndata")
  obj2 <- to_se(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "se"))
})

test_that("to_se works from loom", {
  skip_if_no_loomr()
  library(Seurat)
  obj <- example_obj("loom")
  obj2 <- to_se(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "se"))
})

test_that("to_se works from CDS", {
  skip_if_no_monocle3()
  obj <- example_obj("cds")
  obj2 <- to_se(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "se"))
})
