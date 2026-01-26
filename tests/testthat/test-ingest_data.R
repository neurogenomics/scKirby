test_that("ingest_data defaults to SCE output", {
  obj <- example_obj("seurat")
  obj2 <- ingest_data(obj, verbose = FALSE)
  testthat::expect_true(methods::is(obj2, "SingleCellExperiment"))
})

test_that("ingest_data works with output_type Seurat", {
  obj <- example_obj("sce")
  obj2 <- ingest_data(obj, output_class = "Seurat", verbose = FALSE)
  testthat::expect_true(is_class(obj2, "seurat"))
})

test_that("ingest_data works from list", {
  obj <- example_obj("list")
  obj2 <- ingest_data(obj, verbose = FALSE)
  testthat::expect_true(methods::is(obj2, "SingleCellExperiment"))
})

test_that("ingest_data works from matrix", {
  obj <- example_obj("matrix")
  obj2 <- ingest_data(obj, verbose = FALSE)
  testthat::expect_true(methods::is(obj2, "SingleCellExperiment"))
})

test_that("ingest_data works from CTD", {
  obj <- example_obj("ctd")
  obj2 <- ingest_data(obj, verbose = FALSE)
  testthat::expect_true(methods::is(obj2, "SingleCellExperiment"))
})

test_that("ingest_data works from h5Seurat", {
  skip_if_no_seuratdisk()
  obj <- example_obj("h5seurat")
  obj2 <- ingest_data(obj, verbose = FALSE)
  testthat::expect_true(methods::is(obj2, "SingleCellExperiment"))
})

test_that("ingest_data with return_save_path returns list", {
  obj <- example_obj("seurat")
  save_path <- tempfile(fileext = ".rds")
  result <- ingest_data(obj,
                        save_path = save_path,
                        return_save_path = TRUE,
                        verbose = FALSE)
  testthat::expect_true(is.list(result))
  testthat::expect_true("save_path" %in% names(result))
  testthat::expect_true("obj" %in% names(result))
})
