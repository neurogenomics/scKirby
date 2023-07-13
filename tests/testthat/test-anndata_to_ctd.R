test_that("anndata_to_ctd works", {

  obj <- example_obj("anndata")
  testthat::expect_error(
    anndata_to_ctd(obj)
  )
  obj2 <- anndata_to_ctd(obj,
                         annotLevels = list(groups=NULL),
                         chunk_size = 50)
  testthat::expect_true(is_class(obj2,"ewce"))

  obj3 <- anndata_to_ctd(obj,
                         annotLevels = list(groups=obj$obs$groups))
  testthat::expect_true(is_class(obj3,"ewce"))
})
