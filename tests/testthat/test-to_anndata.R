test_that("to_anndata works from list", {
  skip_if_no_anndata()
  obj <- example_obj("list")
  obj2 <- to_anndata(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "anndata"))
})

test_that("to_anndata works from Seurat", {
  skip_if_no_anndata()
  obj <- example_obj("seurat")
  obj2 <- to_anndata(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "anndata"))
})

test_that("to_anndata works from SE", {
  skip_if_no_anndata()
  obj <- example_obj("se")
  obj2 <- to_anndata(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "anndata"))
})

test_that("to_anndata returns anndata object as-is", {
  skip_if_no_anndata()
  obj <- example_obj("anndata")
  obj2 <- to_anndata(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "anndata"))
})

test_that("to_anndata works from CDS", {
  skip_if_no_anndata()
  skip_if_no_monocle3()
  obj <- example_obj("cds")
  obj2 <- to_anndata(obj, verbose = FALSE)
  testthat::expect_true(is_class(obj2, "anndata"))
})
