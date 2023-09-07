test_that("calc_cor works", {

  obj <- example_obj("seurat")
  Xcor <- calc_cor(obj)
  testthat::expect_true(methods::is(Xcor,"matrix"))
})
