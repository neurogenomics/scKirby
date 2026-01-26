# Test helper functions for conditional skipping

skip_if_no_python <- function() {
  if (!reticulate::py_available(initialize = FALSE)) {
    testthat::skip("Python not available")
  }
}

skip_if_no_anndata <- function() {
  skip_if_no_python()
  if (!reticulate::py_module_available("anndata")) {
    testthat::skip("anndata module not available")
  }
}

skip_if_no_seuratdisk <- function() {
  if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
    testthat::skip("SeuratDisk not installed")
  }
}

skip_if_no_seuratwrappers <- function() {
  if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
    testthat::skip("SeuratWrappers not installed")
  }
}

skip_if_no_monocle3 <- function() {
  if (!requireNamespace("monocle3", quietly = TRUE)) {
    testthat::skip("monocle3 not installed")
  }
}

skip_if_no_loomr <- function() {
  if (!requireNamespace("loomR", quietly = TRUE)) {
    testthat::skip("loomR not installed")
  }
}
