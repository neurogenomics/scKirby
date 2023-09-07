test_that("set_graph works", {

  obj <- example_obj("seurat")
  g <- get_cor(obj)
  key <- "cor_graph"
  obj2 <- set_graph(obj=obj, g=g, key=key)

  g <- get_graphs(obj2, keys = key)
  testthat::expect_true(methods::is(g$cor_graph,"Graph"))
  testthat::expect_true(is_class(g$cor_graph,"matrix"))
})
