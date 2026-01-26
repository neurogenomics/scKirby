test_that("save_data works with RDS format", {
  obj <- example_obj("seurat")
  save_path <- tempfile(fileext = ".rds")
  result <- save_data(obj, filetype = "rds", save_path = save_path, verbose = FALSE)
  testthat::expect_true(file.exists(result))
  # Verify file can be read back
  obj2 <- readRDS(result)
  testthat::expect_true(is_class(obj2, "seurat"))
})

test_that("save_data works with RDA format", {
  obj <- example_obj("seurat")
  save_path <- tempfile(fileext = ".rda")
  result <- save_data(obj, filetype = "rda", save_path = save_path, verbose = FALSE)
  testthat::expect_true(file.exists(result))
})

test_that("save_data works with CDS object", {
  skip_if_no_monocle3()
  obj <- example_obj("cds")
  save_path <- tempfile(fileext = ".rds")
  result <- save_data(obj, filetype = "rds", save_path = save_path, verbose = FALSE)
  testthat::expect_true(file.exists(result))
})

test_that("save_data works with SCE object", {
  obj <- example_obj("sce")
  save_path <- tempfile(fileext = ".rds")
  result <- save_data(obj, filetype = "rds", save_path = save_path, verbose = FALSE)
  testthat::expect_true(file.exists(result))
})

test_that("save_data works with h5 format for SCE", {
  obj <- example_obj("sce")
  save_path <- file.path(tempdir(), "test_hdf5se")
  result <- save_data(obj, filetype = "h5", save_path = save_path, verbose = FALSE)
  testthat::expect_true(dir.exists(result) || file.exists(result))
})

test_that("save_data works with h5seurat format", {
  skip_if_no_seuratdisk()
  obj <- example_obj("seurat")
  save_path <- tempfile(fileext = ".h5seurat")
  result <- save_data(obj, filetype = "h5seurat", save_path = save_path, verbose = FALSE)
  testthat::expect_true(file.exists(result))
})

test_that("save_data works with anndata format", {
  skip_if_no_anndata()
  obj <- example_obj("sce")
  save_path <- tempfile(fileext = ".h5ad")
  result <- save_data(obj, filetype = "anndata", save_path = save_path, verbose = FALSE)
  testthat::expect_true(file.exists(result))
})
