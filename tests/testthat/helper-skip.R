# Test helper functions for conditional skipping

skip_if_no_python <- function() {
  if (!reticulate::py_available(initialize = FALSE)) {
    testthat::skip("Python not available")
  }
}

skip_if_no_anndata <- function() {
  # Check R anndata package is available with working basilisk
  if (!requireNamespace("anndata", quietly = TRUE)) {
    testthat::skip("anndata R package not installed")
  }
  # Check basilisk.utils compatibility (required by anndata)
  if (!requireNamespace("basilisk.utils", quietly = TRUE) ||
      !"activateEnvironment" %in% getNamespaceExports("basilisk.utils")) {
    testthat::skip("basilisk.utils activateEnvironment not available")
  }
  # Also check Python anndata module
  skip_if_no_python()
  if (!reticulate::py_module_available("anndata")) {
    testthat::skip("anndata Python module not available")
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
