test_that("list_to_anndata works", {

  #### From obj ####
  obj <- example_obj("list")
  obj2 <- list_to_anndata(obj)
  testthat::expect_true(is_class(obj2,"anndata"))

  #### From path names ####
  obj3 <- example_obj("list_paths")
  obj4 <- list_to_anndata(obj = obj3)
  testthat::expect_true(is_class(obj4,"anndata"))
})
