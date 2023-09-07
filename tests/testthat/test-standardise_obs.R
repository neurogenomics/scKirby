test_that("standardise_obs works", {

  level_cols <- list(level1="")

  obj <- example_obj("anndata")
  obj2 <- standardise_obs(obj = obj,
                           level_cols = level_cols)
})
