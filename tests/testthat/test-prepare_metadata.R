test_that("prepare_metadata works", {

  level_cols <- list(level1="")

  obj <- example_obj("anndata")
  obj2 <- prepare_metadata(obj = obj,
                           level_cols = level_cols)
})
