test_that("get_cor works", {

  obj <- example_obj("seurat")
  obj2 <- get_cor(obj = obj,
                  keys = "pca",
                  return_obj = TRUE)
  testthat::expect_true(is_class(obj2,"seurat"))
  Xcor <- get_cor(obj = obj,
                  keys = "pca",
                  return_obj = FALSE)
  testthat::expect_true(methods::is(Xcor,"Graph"))
  testthat::expect_true(is_class(Xcor,"matrix"))

  Xcor2 <- get_cor(obj = obj,
                   keys = "pca",
                   return_obj = FALSE,
                   as_graph = FALSE)
  testthat::expect_true(methods::is(Xcor2,"matrix"))
  testthat::expect_true(is_class(Xcor2,"matrix"))
})
