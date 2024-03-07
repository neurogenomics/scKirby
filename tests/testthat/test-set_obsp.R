test_that("set_obsp works", {

  obj <- example_obj("seurat")
  obsp <- get_cor(obj)
  key <- "cor_graph"
  obj2 <- set_obsp(obj=obj, obsp=obsp, key=key)

  obsp <- get_obsp(obj2, keys = key)
  testthat::expect_true(methods::is(obsp$cor_graph,"Graph"))
  testthat::expect_true(is_class(obsp$cor_graph,"matrix"))
})
