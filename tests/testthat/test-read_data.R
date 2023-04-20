test_that("read_data works", {

  #### RDS: seurat ####
  obj <- example_obj("seurat")
  tmp <- tempfile("seurat")
  saveRDS(obj,tmp)
  obj2 <- read_data(tmp)
  testthat::expect_true(is_class(obj2,"seurat"))
  #### RDS: cds ####
  obj <- example_obj("cds")
  tmp <- tempfile("cds")
  saveRDS(obj,tmp)
  obj2 <- read_data(tmp)
  testthat::expect_true(is_class(obj2,"cds"))
  #### RDS: sce ####
  obj <- example_obj("sce")
  tmp <- tempfile("sce")
  saveRDS(obj,tmp)
  obj2 <- read_data(tmp)
  testthat::expect_true(is_class(obj2,"se"))
  #### anndata ####
  tmp <- tempfile(fileext = "anndata.h5ad")
  obj <- example_obj("anndata", save_path = tmp)
  obj2 <- read_data(tmp)
  testthat::expect_true(is_class(obj2,"anndata"))
  #### loom ####
  library(Seurat)
  obj <- example_obj("loom")
  obj2 <- read_data(obj$filename)
  testthat::expect_true(is_class(obj2,"loom"))
  #### Matrix/data.table ####
  obj <- example_obj("data.table")
  tmp <- tempfile(fileext = ".tsv.gz")
  data.table::fwrite(obj,file = tmp,
                     sep = "\t")
  obj2 <- read_data(path = tmp)
  testthat::expect_true(is_class(obj2,"matrix"))
})
