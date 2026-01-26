test_that("get_obs works for Seurat", {
  obj <- example_obj("seurat")
  obs <- get_obs(obj, verbose = FALSE)
  testthat::expect_true(is.data.frame(obs) || methods::is(obs, "DFrame"))
})

test_that("get_obs works for SE", {
  obj <- example_obj("se")
  obs <- get_obs(obj, verbose = FALSE)
  testthat::expect_true(is.data.frame(obs) || methods::is(obs, "DFrame"))
})

test_that("get_obs works for SCE", {
  obj <- example_obj("sce")
  obs <- get_obs(obj, verbose = FALSE)
  testthat::expect_true(is.data.frame(obs) || methods::is(obs, "DFrame"))
})

test_that("get_obs works for matrix", {
  obj <- example_obj("matrix")
  obs <- get_obs(obj, verbose = FALSE)
  testthat::expect_true(is.data.frame(obs))
})

test_that("get_obs works for list", {
  obj <- example_obj("list")
  obs <- get_obs(obj, verbose = FALSE)
  testthat::expect_true(is.data.frame(obs) || methods::is(obs, "DFrame"))
})

test_that("get_obs works for h5Seurat", {
  skip_if_no_seuratdisk()
  obj <- example_obj("h5seurat")
  obs <- get_obs(obj, verbose = FALSE)
  testthat::expect_true(is.data.frame(obs) || methods::is(obs, "DFrame"))
})

test_that("get_obs works for anndata", {
  skip_if_no_anndata()
  obj <- example_obj("anndata")
  obs <- get_obs(obj, verbose = FALSE)
  testthat::expect_true(is.data.frame(obs) || methods::is(obs, "DFrame"))
})

test_that("get_obs works for loom", {
  skip_if_no_loomr()
  library(Seurat)
  obj <- example_obj("loom")
  obs <- get_obs(obj, verbose = FALSE)
  testthat::expect_true(is.data.frame(obs))
})

test_that("get_obs works for CDS", {
  skip_if_no_monocle3()
  obj <- example_obj("cds")
  obs <- get_obs(obj, verbose = FALSE)
  testthat::expect_true(is.data.frame(obs) || methods::is(obs, "DFrame"))
})
