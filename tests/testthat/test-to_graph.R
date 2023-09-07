test_that("to_graph works", {

  obj <- example_obj("matrix")
  obj2 <- to_graph(obj)
  testthat::expect_true(methods::is(obj2,"Graph"))
})
