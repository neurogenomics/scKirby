test_that("get_x works for Seurat", {
  obj <- example_obj("seurat")
  X <- get_x(obj, verbose = FALSE)
  testthat::expect_true(is.list(X) || is.matrix(X) || methods::is(X, "sparseMatrix"))
})

test_that("get_x works for SE", {
  obj <- example_obj("se")
  X <- get_x(obj, verbose = FALSE)
  testthat::expect_true(is.list(X) || is.matrix(X) || methods::is(X, "sparseMatrix"))
})

test_that("get_x works for SCE", {
  obj <- example_obj("sce")
  X <- get_x(obj, verbose = FALSE)
  testthat::expect_true(is.list(X) || is.matrix(X) || methods::is(X, "sparseMatrix"))
})

test_that("get_x works for matrix", {
  obj <- example_obj("matrix")
  X <- get_x(obj, verbose = FALSE)
  testthat::expect_true(is.matrix(X) || methods::is(X, "sparseMatrix"))
})

test_that("get_x works for list", {
  obj <- example_obj("list")
  X <- get_x(obj, verbose = FALSE)
  testthat::expect_true(is.list(X) || is.matrix(X) || methods::is(X, "sparseMatrix"))
})

test_that("get_x with n parameter returns limited elements", {
  obj <- example_obj("seurat")
  X <- get_x(obj, n = 1, verbose = FALSE)
  testthat::expect_true(length(X) <= 1 || !is.list(X))
})

test_that("get_x with transpose=TRUE transposes data", {
  obj <- example_obj("seurat")
  X1 <- get_x(obj, n = 1, verbose = FALSE)[[1]]
  X2 <- get_x(obj, n = 1, transpose = TRUE, verbose = FALSE)[[1]]
  # Dimensions should be swapped
  testthat::expect_equal(nrow(X1), ncol(X2))
  testthat::expect_equal(ncol(X1), nrow(X2))
})

test_that("get_x works for h5Seurat", {
  skip_if_no_seuratdisk()
  obj <- example_obj("h5seurat")
  X <- get_x(obj, verbose = FALSE)
  testthat::expect_true(is.list(X) || is.matrix(X) || methods::is(X, "sparseMatrix"))
})

test_that("get_x works for anndata", {
  skip_if_no_anndata()
  obj <- example_obj("anndata")
  X <- get_x(obj, verbose = FALSE)
  testthat::expect_true(is.list(X) || is.matrix(X) || methods::is(X, "sparseMatrix"))
})

test_that("get_x works for loom", {
  skip_if_no_loomr()
  library(Seurat)
  obj <- example_obj("loom")
  X <- get_x(obj, verbose = FALSE)
  testthat::expect_true(is.list(X) || is.matrix(X) || methods::is(X, "sparseMatrix"))
})

test_that("get_x works for CDS", {
  skip_if_no_monocle3()
  obj <- example_obj("cds")
  X <- get_x(obj, verbose = FALSE)
  testthat::expect_true(is.list(X) || is.matrix(X) || methods::is(X, "sparseMatrix"))
})
